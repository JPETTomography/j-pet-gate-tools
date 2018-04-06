#ifndef GateAttenuationCorrectionActorMESSENGER_CC
#define GateAttenuationCorrectionActorMESSENGER_CC

#include "GateAttenuationCorrectionActorMessenger.hh"
#include "GateAttenuationCorrectionActor.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//-----------------------------------------------------------------------------
GateAttenuationCorrectionActorMessenger::GateAttenuationCorrectionActorMessenger(GateAttenuationCorrectionActor* sensor):GateActorMessenger(sensor),pGateAttenuationCorrectionActor(sensor)
{
	BuildCommands(baseName+sensor->GetObjectName());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
GateAttenuationCorrectionActorMessenger::~GateAttenuationCorrectionActorMessenger()
{
	delete pSetNumProjectionBinsCmd;
	delete pSetNumAngleBinsCmd;
	delete pSetNumAxialSinogramsCmd;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateAttenuationCorrectionActorMessenger::BuildCommands(G4String base)
{
	G4String guidance;
	G4String bb;
	
	bb = base+"/setNumProjectionBins";
	pSetNumProjectionBinsCmd = new G4UIcmdWithAnInteger(bb,this);
	guidance = "Set the number of projection bins in the sinogram.";
	pSetNumProjectionBinsCmd->SetGuidance(guidance);
	
	bb = base+"/setNumAngleBins";
	pSetNumAngleBinsCmd = new G4UIcmdWithAnInteger(bb,this);
	guidance = "Set the number of angle bins in the sinogram.";
	pSetNumAngleBinsCmd->SetGuidance(guidance);
	
	bb = base+"/setNumAxialSinograms";
	pSetNumAxialSinogramsCmd = new G4UIcmdWithAnInteger(bb,this);
	guidance = "Set the number of axial sinograms (GateAttenuationCorrectionActor only procudes 2D slices).";
	pSetNumAxialSinogramsCmd->SetGuidance(guidance);
	
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateAttenuationCorrectionActorMessenger::SetNewValue(G4UIcommand* command, G4String param)
{
	if(command == pSetNumProjectionBinsCmd) pGateAttenuationCorrectionActor->SetNumProjectionBins(pSetNumProjectionBinsCmd->GetNewIntValue(param));
	if(command == pSetNumAngleBinsCmd) pGateAttenuationCorrectionActor->SetNumAngleBins(pSetNumAngleBinsCmd->GetNewIntValue(param));
	if(command == pSetNumAxialSinogramsCmd) pGateAttenuationCorrectionActor->SetNumAxialSinograms(pSetNumAxialSinogramsCmd->GetNewIntValue(param));
	
	
	GateActorMessenger::SetNewValue(command ,param );
}
//-----------------------------------------------------------------------------

#endif /* end #define GateAttenuationCorrectionActorMESSENGER_CC */
