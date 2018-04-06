/*
 \class  GateAttenuationCorrectionActorMessenger
 */

#ifndef GateAttenuationCorrectionActorMESSENGER_HH
#define GateAttenuationCorrectionActorMESSENGER_HH

#include "globals.hh"
#include "GateActorMessenger.hh"

class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class GateAttenuationCorrectionActor;
class GateAttenuationCorrectionActorMessenger : public GateActorMessenger
{
public:
	GateAttenuationCorrectionActorMessenger(GateAttenuationCorrectionActor* sensor);
	virtual ~GateAttenuationCorrectionActorMessenger();
	
	void BuildCommands(G4String base);
	void SetNewValue(G4UIcommand*, G4String);
	
protected:
	GateAttenuationCorrectionActor * pGateAttenuationCorrectionActor;
	
	G4UIcmdWithAnInteger * pSetNumProjectionBinsCmd;
	G4UIcmdWithAnInteger * pSetNumAngleBinsCmd;
	G4UIcmdWithAnInteger * pSetNumAxialSinogramsCmd;
};

#endif /* end #define GateAttenuationCorrectionActorMESSENGER_HH*/
