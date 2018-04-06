/*!
 \class GateAttenuationCorrectionActor
 \author gbindsei@uwo.ca
 */

#ifndef GATEAttenuationCorrectionActor_HH
#define GATEAttenuationCorrectionActor_HH

#include "GateVActor.hh"
#include "GateAttenuationCorrectionActorMessenger.hh"
#include "GateActorManager.hh"
#include "G4EmCalculator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"

// struct to store two G4ThreeVectors
struct ParticleInitialConditions
{
	G4ThreeVector position;
	G4ThreeVector direction;
};

//-----------------------------------------------------------------------------
/// \brief Actor computing attenuation correction factors for LORs of a ECAT-like PET system
class GateAttenuationCorrectionActor : public GateVActor
{

private:

	G4double PI; // Definition of the number PI
    
	size_t num_of_rings; // The number of rings
	size_t num_of_crystals_in_ring; // The number of crystals in a single ring	
	size_t num_projection_bins; // The number of projection bins in the sinogram
	size_t num_angle_bins; // The number of angle bins in the sinogram
	size_t num_axial_sinograms; // The number of axial sinograms that make up the 3D sinogram
	
	G4bool bHasEnteredFOV; // True if the geantino has entered the FOV (world volume) or false if it hasn't yet left the first crystal block
    G4double mGlobalACScaleFactor; // Factor to scale each linear attenuation coefficent (optional).
	G4double mAirAttenuationValue; // Linear attenuation coefficient of the world default air.
	
    size_t	theRingNum1, theRingNum2, theCrystalNum1, theCrystalNum2; // Variables to store the current crystal identifiers
	
	G4String mCrystalVolumeName; // Stores the name of the crystal volume (to determine if the step point is within a crystal)
	
	// Temporarily stores the running sum of the attenuation coefficient times the step distance.
	// At the end of an event, the attenuation correction factor, ACF = exp( - mMuTimesDeistanceRunningSum );
	G4double mMuTimesDistanceRunningSum;
	
	G4ThreeVector** mCrystalPositionArray; // Declare a pointer to a pointer array to the array that will store the crystal positions (size is determined at run time)
	
	G4double*** mAttenuationCorrectionFactorSinogram; // Declare a triple pointer to an array that will store the attenuation correction factors (size set at run time)
	
	void ComputeCrystalPositions();
	void ComputeAttenuationCoefficients();
	
	ParticleInitialConditions CalculateParticleInitialConditions(size_t aRingNumInitial, size_t aCrystalNumInitial, size_t aRingNumFinal, size_t aCrystalNumFinal);
	
	G4double* mMaterialAttenuationCoefficientArray; // Declare a pointer to the array of attenuation coefficients for each material
	
public: 
	
	virtual ~GateAttenuationCorrectionActor();
    
	//-----------------------------------------------------------------------------
	// This macro initialize the CreatePrototype and CreateInstance
	FCT_FOR_AUTO_CREATOR_ACTOR(GateAttenuationCorrectionActor)
	
	//-----------------------------------------------------------------------------
	// Constructs the sensor
	virtual void Construct();
	
	//-----------------------------------------------------------------------------
	// Callbacks
    virtual void BeginOfRunAction(const G4Run*);
	virtual void BeginOfEventAction(const G4Event*);
	virtual void EndOfEventAction(const G4Event*);
	//    virtual void PreUserTrackingAction(const GateVVolume *, const G4Track*);
    virtual void UserSteppingAction(const GateVVolume *, const G4Step*); 
	
	// other methods
	void SetEnergy (double E) {mEnergy=E;}
	void SetParticleName (G4String Name) {mPartName=Name;}
	void SetAttenFactor (size_t bin_number, G4double theACF);
	void SetNumProjectionBins (G4int aNumProjectionBins) {num_projection_bins = aNumProjectionBins;}
	void SetNumAngleBins (G4int aNumAngleBins) {num_angle_bins = aNumAngleBins;}
	void SetNumAxialSinograms (G4int aNumAxialSinograms) {num_axial_sinograms = aNumAxialSinograms;}
	
	// methods to return the global position of a crystal given the ring number and crystal ID
	// ring numbers start from zero at the ring farthest from the bed
	// looking from the bed side, the zeroth crystal is the top most (slightly to the right)
	// the crystal ID's go from 0-319 and the ring numbers go from 0-79
	inline G4ThreeVector GetGlobalCrystalPosition(size_t ring_number, size_t crystal_id)
	{
		return mCrystalPositionArray[ring_number][crystal_id];
	}
	//-----------------------------------------------------------------------------
	/// Saves the data collected to the file
	virtual void SaveData();
	virtual void ResetData();
	
protected:
	
	double mEnergy;
	G4String mPartName;
	
	G4EmCalculator * emcalc;
	GateAttenuationCorrectionActor(G4String name, G4int depth=0); 
	GateAttenuationCorrectionActorMessenger * pActorMessenger;
	//  GateActorMessenger * pActor;
};

MAKE_AUTO_CREATOR_ACTOR(AttenuationCorrectionActor,GateAttenuationCorrectionActor)


#endif /* end #define GateAttenuationCorrectionActor_HH */

