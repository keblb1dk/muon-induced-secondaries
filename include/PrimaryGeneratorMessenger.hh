#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UicmdWithADouble;

class PrimaryGeneratorMessenger: public G4UImessenger
{
public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
    virtual ~PrimaryGeneratorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
private:
    PrimaryGeneratorAction* fPrimaryAction;
    
    G4UIdirectory*          fPrimaryDir;
    G4UIcmdWithABool*       fUseMeanSamplingCmd;
    G4UIcmdWithAString*     fAngularFileCmd;
    G4UIcmdWithAString*     fMeanEnergyFileCmd;
    G4UIcmdWithAString*     fCDFEnergyFileCmd;
    G4UIcmdWithADoubleAndUnit*     fEnergyThresholdCmd;
};

#endif
