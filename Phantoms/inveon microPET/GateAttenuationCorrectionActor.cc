/*
 \brief Class GateAttenuationCorrectionActor : 
 \brief 
 */

#ifndef GateAttenuationCorrectionActor_CC
#define GateAttenuationCorrectionActor_CC

#include "GateAttenuationCorrectionActor.hh"
#include "GateMiscFunctions.hh"
#include "G4Event.hh"
#include "G4MaterialTable.hh"
#include "G4ParticleTable.hh"
#include "G4UnitsTable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "GateVolumeID.hh"
#include "GateVVolume.hh"
#include "GateVSystem.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4StepPoint.hh"
#include "GateSystemListManager.hh"
#include "GateSystemComponent.hh"
#include "GateLinearRepeater.hh"
#include "GateBox.hh"

//-----------------------------------------------------------------------------
/// Constructors (Prototype)
GateAttenuationCorrectionActor::GateAttenuationCorrectionActor(G4String name, G4int depth):
GateVActor(name,depth)
{
	GateDebugMessageInc("Actor",4,"GateAttenuationCorrectionActor() -- begin"<<G4endl);  
	//SetTypeName("EmCalculatorActor");
	//  pActor = new GateActorMessenger(this);
	ResetData();
	GateDebugMessageDec("Actor",4,"GateAttenuationCorrectionActor() -- end"<<G4endl);
	
	mEnergy = 511.*keV;
	mPartName = "gamma";
	
	num_projection_bins = 128;
	num_angle_bins = 160;
	num_axial_sinograms = 159;
	
    mGlobalACScaleFactor = 1/1.028; // This scale factor reduces the ACF values to remove the effect of small-angle Compton scattering.
	mAirAttenuationValue = 0.;
    
    PI = atan(1.0)*4;

	pActorMessenger = new GateAttenuationCorrectionActorMessenger(this);
	emcalc = new G4EmCalculator;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/// Destructor 
GateAttenuationCorrectionActor::~GateAttenuationCorrectionActor() 
{
	delete pActorMessenger;
	// Deallocate mAttenuationCorrectionFactorSinogram
	for (size_t i = 0; i < num_projection_bins; ++i)
	{
		for (size_t j = 0; j < num_angle_bins; ++j)
		{
			delete [] mAttenuationCorrectionFactorSinogram[i][j];
		}
		delete [] mAttenuationCorrectionFactorSinogram[i];
	}
	delete [] mAttenuationCorrectionFactorSinogram;
	
	// Deallocate mCrystalPositionArray;
	for (size_t i = 0; i < num_of_rings; ++i)
	{
		delete [] mCrystalPositionArray[i];
	}
	delete [] mCrystalPositionArray;
	
	delete [] mMaterialAttenuationCoefficientArray;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/// Construct
void GateAttenuationCorrectionActor::Construct()
{
	GateVActor::Construct();
	//  Callbacks
	EnableBeginOfRunAction(true);
	EnableBeginOfEventAction(true);
	EnableEndOfEventAction(true);
	EnableUserSteppingAction(true);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// SetAttenFactor
void GateAttenuationCorrectionActor::SetAttenFactor(size_t bin_number, G4double theACF)
{
	// This method sets the bin number in the Attenuation Correction Factor sinogram to the value theACF
	G4int projection_index = bin_number%num_projection_bins;
	G4int angle_index = (bin_number/num_projection_bins)%num_angle_bins;
	G4int axial_index = bin_number/(num_projection_bins*num_angle_bins);	
		
	mAttenuationCorrectionFactorSinogram[projection_index][angle_index][axial_index] = theACF;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// CalculateParticleInitialConditions
ParticleInitialConditions GateAttenuationCorrectionActor::CalculateParticleInitialConditions(size_t aRingNumInitial, size_t aCrystalNumInitial, size_t aRingNumFinal, size_t aCrystalNumFinal)
{
	// This function returns the position and direction for a particle along the line of response between two crystals
	// Arguments: ring and crystal indices for the initial position, ring and crystal indicies for the final position
	// Output: struct containing G4ThreeVectors for position (initial) and direction (towards the final position)
	
	ParticleInitialConditions theInitialConditions;
	
	theInitialConditions.position = mCrystalPositionArray[aRingNumInitial][aCrystalNumInitial];
		
	theInitialConditions.direction = (mCrystalPositionArray[aRingNumFinal][aCrystalNumFinal] - mCrystalPositionArray[aRingNumInitial][aCrystalNumInitial]);
	
	return theInitialConditions;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ComputeCrystalPositions
void GateAttenuationCorrectionActor::ComputeCrystalPositions()
{
	// This method fills mCrystalPositionArray with G4ThreeVector elements containing
	// the positions of the centres of all the crystal elements in the detector.
	// 
	// ring numbers start from zero at the ring farthest from the bed
	// looking from the bed side, the zeroth crystal is the top most (slightly to the right)
	// the crystal ID's go from 0-319 and the ring numbers go from 0-79
    // The zeroth crystal should be the leftmost crystal in the top right block, looking from bed side.
	
	// Get a pointer to the detector system
	GateSystemListManager* theSystemList = GateSystemListManager::GetInstance();
	GateVSystem* theSystem = theSystemList->GetSystem(0);
	
	// Get pointers to the Base, Block and Crystal components.
	GateSystemComponent* theBaseComponent = theSystem->GetBaseComponent();
	GateSystemComponent* theBlockComponent = theBaseComponent->GetChildComponent(0);
	GateSystemComponent* theCrystalComponent = theBlockComponent->GetChildComponent(0);
	
	G4ThreeVector aPositionCrystal, aPositionBlock, aPositionSystem;
	G4RotationMatrix* aRotationBlockPtr;
	G4RotationMatrix aRotationBlock;
	G4RotationMatrix aIdentityMatrix;
	
	mCrystalVolumeName = theCrystalComponent->GetPhysicalVolume()->GetName();
	
	// Get the geometry of the ECAT-like system
	GateArrayComponent* theCrystalArrayComponent = theSystem->FindArrayComponent("crystal");
	
	size_t B_z = theBlockComponent->GetLinearRepeatNumber(); // Axial number of blocks (4 for Inveon)
    size_t B_phi = theBlockComponent->GetAngularRepeatNumber(); // Azimuthal number of blocks (16 for Inveon)
    size_t C_z = theCrystalArrayComponent->GetRepeatNumber(2); // Axial number of crystals in z direction (20 for Inveon)
    size_t C_phi;
	// Check to see which second axis (x or y) the user chose to apply the repeater on.
	if (theCrystalArrayComponent->GetRepeatVector().x() == 0)
	{
		C_phi = theCrystalArrayComponent->GetRepeatNumber(1); // Tangential number of crystals in y-direction (20 for Inv)
	}
	else
	{
		C_phi = theCrystalArrayComponent->GetRepeatNumber(0); // Tangential number of crystals in x-direction (20 for Inv)
	}
	
	num_of_rings = B_z*C_z;
	num_of_crystals_in_ring = B_phi*C_phi;
	
	// Initialize the dynamic array mCrystalPositionArray
	// This code is more complicated than a usual dynamic array initialization because
	// mCrystalPositionArray should have two dimensions: mCrystalPositionArray[num_of_rings][num_of_crystals_in_ring]
	mCrystalPositionArray = new G4ThreeVector* [num_of_rings]; // initialize ring index
	for (size_t i = 0; i < num_of_rings; ++i)
	{
		mCrystalPositionArray[i] = new G4ThreeVector [num_of_crystals_in_ring]; // initialize crystal index
	}
	
	// Get the translation of the overall system
	aPositionSystem = theBaseComponent->GetCurrentTranslation(0);
	
	size_t blockID, crystalID;
	
	// Loop over the rings and the number of crystals in each ring to fill the mCrystalPositionArray array
	for (size_t ringindex = 0; ringindex < num_of_rings; ringindex++)
	{
		for (size_t crystalindex = 0; crystalindex < num_of_crystals_in_ring; crystalindex++)
		{
			// Compute the GATE specific crystal and block identifiers
			crystalID = crystalindex%C_phi + C_phi*(ringindex%C_z);
			blockID = B_z*floor(crystalindex/C_phi) + floor(ringindex/C_z);			
			
			aPositionCrystal = theCrystalComponent->GetCurrentTranslation(crystalID);
			aPositionBlock = theBlockComponent->GetCurrentTranslation(blockID);
			aRotationBlockPtr = theBlockComponent->GetCurrentRotation(blockID);
			// If no rotation pointer, then set aRotationBlock to point to an identity matrix
			if (!aRotationBlockPtr) { aRotationBlock = aIdentityMatrix; } else { aRotationBlock = *aRotationBlockPtr; }
			
			// Now reverse the sign of the angle of the rotation matrix (seems to be necessary at least for my crystal definition)
			aRotationBlock.setDelta( - aRotationBlock.getDelta() );
			
            // Note: I shift crystalindex over by 10 to be consistent with the Inveon convention that the 0th crystal bin is at the very top of the scanner (need to shift 10 cyrstals over)            
			// I think that the -(C_phi/2) term does not have to be there if the detectors are rotated properly.
            mCrystalPositionArray[ringindex][crystalindex] = aPositionSystem + aPositionBlock + aRotationBlock*aPositionCrystal;
		}
	}
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Compute the attenuation coefficients material look up array
void GateAttenuationCorrectionActor::ComputeAttenuationCoefficients()
{
	// Get the material table
	const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
	
	// Compute the number of materials in the simulation
	size_t theNumberOfMaterials = G4Material::GetNumberOfMaterials();
	
	//  Initialize the attenuation coefficients material look up array
	mMaterialAttenuationCoefficientArray = new G4double[theNumberOfMaterials];
	
	// Declare temporary variables to store the cross sections and material name.
	G4String aMaterialName;
	G4double xCompton, xPhotoElec, xTotal;
	G4double cut = DBL_MAX;
    
    G4cout << G4endl;
    G4cout << "List of Materials and Attenuation Coefficients for 0.511 MeV Photons:" << G4endl;
    G4cout << "[Material Name] [mu_Compton 1/mm] [mu_PhotoElec 1/mm] [mu_Total 1/mm]" << G4endl;
    
	// Loop over all the materials in the simulation and store the total linear attenuation coefficient for each material
	for(size_t k=0;k<theNumberOfMaterials;k++)
    {
		aMaterialName = (*theMaterialTable)[k]->GetName();
		// ComputeCrossSectionPerVolume returns the linear attenuation coefficient in default units 1/mm
		xCompton = emcalc->ComputeCrossSectionPerVolume(mEnergy, mPartName, "Compton", aMaterialName, cut);
		xPhotoElec = emcalc->ComputeCrossSectionPerVolume(mEnergy, mPartName, "PhotoElectric", aMaterialName, cut);
		xTotal = xCompton + xPhotoElec;
		mMaterialAttenuationCoefficientArray[k] = xTotal;	
		if (aMaterialName == "worldDefaultAir") { mAirAttenuationValue = xTotal; }
        
        // Print out List of Materials and Attenuation Coefficients for 0.511 MeV Photons
        G4cout << aMaterialName << "\t" << xCompton << "\t" << xPhotoElec <<"\t"<< xTotal << G4endl;
    }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Callback Begin of Run
void GateAttenuationCorrectionActor::BeginOfRunAction(const G4Run*r)
{
	G4cout << G4endl;
	G4cout << "Gate Attenuation Correction Actor Started" << G4endl;
	G4cout << G4endl;
	G4cout << "Number of Projection Bins: " << num_projection_bins << G4endl;
	G4cout << "Number of Angle Bins: " << num_angle_bins << G4endl;
	G4cout << "NUmber of Axial Sinograms: " << num_axial_sinograms << G4endl;
	G4cout << G4endl;
	
	// Initialize and fill the array of global crystal positions: mCrystalPositionArray[ring][crystal].
	ComputeCrystalPositions();
	
	// Compute the cross sections (attenuation coefficients) for each material in the simulation and store in a lookup table
	ComputeAttenuationCoefficients();
	
	// Initialize the dynamic array of the attenuation correction factor sinogram
	// mAttenuationCorrectionFactorSinogram should 3 dimensions: mAtten...Sinogram[num_of_projection_bins][num_of_angle_bins][num_of_axial_sinograms]
	
	mAttenuationCorrectionFactorSinogram = new G4double** [num_projection_bins]; // Initialize projection bins
	for (size_t i = 0; i < num_projection_bins; ++i)
	{
		mAttenuationCorrectionFactorSinogram[i] = new G4double* [num_angle_bins]; // Initialize angle bins
		for (size_t j = 0; j < num_angle_bins; ++j)
		{
			mAttenuationCorrectionFactorSinogram[i][j] = new G4double [num_axial_sinograms]; // Initialize axial sinograms
			for (size_t k = 0; k < num_axial_sinograms; ++k)
			{
				mAttenuationCorrectionFactorSinogram[i][j][k] = 0.; // Set all ACF's to zero.
			}
		}
	}
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Callback Begin Event
void GateAttenuationCorrectionActor::BeginOfEventAction(const G4Event*e)
 {
	 size_t theEventID = e->GetEventID();
	 
	 // Reset the bHasEnteredFOV flag to be false (for the start of the current event)
	 bHasEnteredFOV = false;
	 
	 // Reset the mMuTimesDistanceRunningSum variable to zero
	 mMuTimesDistanceRunningSum = 0.;
	 
	 // Get a pointer to the PrimaryVertex which contains the position of the particle generation and gives access to the particle
	 G4PrimaryVertex* thePrimaryVertex = e->GetPrimaryVertex();
	 G4PrimaryParticle* thePrimaryParticle = thePrimaryVertex->GetPrimary(0);
	 
	 // Set the particle to be a geantino (non interacting)
	 G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
	 thePrimaryParticle->SetParticleDefinition(theParticleTable->FindParticle("geantino"));
	 
	 // Compute the crystal and ring numbers of the desired LOR for this event
	 size_t initialC1, initialC2, theAngleIteration, theProjectionIteration; // The initial crystals to start and current angle & projection iteration
	 
	 theRingNum1 = theEventID/(2*num_projection_bins*num_angle_bins);
	 theRingNum2 = theRingNum1+(theEventID/(num_projection_bins*num_angle_bins))%2;
	 
	 // The 'fixed angle' way (fill the sinogram across first as you loop through projections)
	 theAngleIteration = (theEventID/num_projection_bins)%num_angle_bins;
	 theProjectionIteration = theEventID%num_projection_bins;
	 
     // Double check that it should be +10 for each. The 64th projection iteration on 0th angle should go through 0, 160 according to siemens ring index.
	 initialC1 = num_of_crystals_in_ring - num_projection_bins/4 + theAngleIteration;
	 initialC2 = num_of_crystals_in_ring/2 + num_projection_bins/4 + theAngleIteration;
	 
	 // Force positive cyrstal numbers
	 theCrystalNum1 = (initialC1 + theProjectionIteration/2 + num_of_crystals_in_ring)%num_of_crystals_in_ring;
	 
	 theCrystalNum2 = (initialC2 - (theProjectionIteration+1)/2 + num_of_crystals_in_ring)%num_of_crystals_in_ring;
	 
	 // Compute the initial position and direction of the particle
	 ParticleInitialConditions theInitialConditions = CalculateParticleInitialConditions(theRingNum1,theCrystalNum1,theRingNum2,theCrystalNum2);
	 thePrimaryVertex->SetPosition(theInitialConditions.position.x(), theInitialConditions.position.y(), theInitialConditions.position.z());
	 thePrimaryParticle->SetMomentum(theInitialConditions.direction.x(), theInitialConditions.direction.y(), theInitialConditions.direction.z());
 }
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Callback End Event
void GateAttenuationCorrectionActor::EndOfEventAction(const G4Event*e)
 {
	 // Double check that the flag bHasEnteredFOV is true (this should always be true at the end of any event or else something is wrong)
	 if (bHasEnteredFOV == false) { G4cerr << "ERROR, an event LOR has not passed through the FOV. Something is wrong." << G4endl; }
	 
	 // Compute the attenuation correction factor (ACF) for the current LOR.
     size_t event_number = e->GetEventID();
          
     // This method sets the bin number in the Attenuation Correction Factor sinogram to the value theACF
     G4int projection_index = event_number%num_projection_bins;
     G4int axial_index = event_number/(num_projection_bins*num_angle_bins);
     G4int angle_index = (event_number/num_projection_bins)%num_angle_bins;
     
     mAttenuationCorrectionFactorSinogram[projection_index][angle_index][axial_index] = exp( mMuTimesDistanceRunningSum * mGlobalACScaleFactor );
	 
     if (event_number % 100000 == 0) { G4cout << "LORs Processed:  " << event_number << G4endl; }
 }
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Callback Step Action
void GateAttenuationCorrectionActor::UserSteppingAction(const GateVVolume * v, const G4Step * step)
{	
	// Get the pointer to the pre step point and get its  material pointer
	G4StepPoint* thePrePoint = step->GetPreStepPoint();
	G4Material* thePreMaterial = thePrePoint->GetMaterial();
	
	// Get the pointer to the post step point
	G4StepPoint* thePostPoint = step->GetPostStepPoint();
	
	// Check if the geantino has entered the FOV (world volume with no mother volume) for the first time and if so, set bHasEnteredFOV to be true
	// Purpose of this code is to ensure a particle is only killed after it has left the original crystal array (it may hit multiple crystals on the way out)
	// And to only add attenuation correction factors 
	if (bHasEnteredFOV == false && thePrePoint->GetPhysicalVolume()->GetMotherLogical() == NULL) { bHasEnteredFOV = true; }
	
	// If the geantino will enter an opposing crystal in the next step, kill the track (still finishes adding the attenuation correction for current step)
	if (bHasEnteredFOV == true && thePostPoint->GetPhysicalVolume()->GetName() == mCrystalVolumeName) { step->GetTrack()->SetTrackStatus( fStopAndKill ); }
	
	// If the geantino is in the FOV and traversing the phantom, accumulate the attenuation amount for the current step.
	// The attenuation amount is stored in mMuTimesDistanceRunningSum and the total ACF = exp (-mMuTimesDistanceRunningSum).
	// Also do not count the attenuation of the world air (since that is accounted for in normalization
	if (bHasEnteredFOV)
	{
		mMuTimesDistanceRunningSum += (mMaterialAttenuationCoefficientArray[thePreMaterial->GetIndex()] - mAirAttenuationValue ) * step->GetStepLength();
	}
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/// Save data
void GateAttenuationCorrectionActor::SaveData()
{
	std::ofstream os;
	os.open(mSaveFilename.data(),std::ios::out | std::ios::binary);
	
	if (!os) {
		GateMessage("Output",1,"Error Writing file: " <<mSaveFilename << G4endl);
	}
	
	G4float temp_ACF;
	for (size_t k = 0; k<num_axial_sinograms; k++)
	{
		for (size_t j = 0; j<num_angle_bins; j++)
		{
			for (size_t i = 0; i<num_projection_bins; i++)
			{
				temp_ACF = static_cast<G4float> (mAttenuationCorrectionFactorSinogram[i][j][k]);
				os.write((char *) &temp_ACF, 4);
			}
		}
	}
	os.close();
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void GateAttenuationCorrectionActor::ResetData() 
{
}
//-----------------------------------------------------------------------------

#endif /* end #define GateAttenuationCorrectionActor_CC */

