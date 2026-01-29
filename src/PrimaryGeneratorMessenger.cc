#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* primaryAction)
 : G4UImessenger(),
   fPrimaryAction(primaryAction),
   fPrimaryDir(0),
   fUseMeanSamplingCmd(0),
   fAngularFileCmd(0),
   fMeanEnergyFileCmd(0),
   fCDFEnergyFileCmd(0),
   fEnergyThresholdCmd(0)
{
    fPrimaryDir = new G4UIdirectory("/primary/");
    fPrimaryDir->SetGuidance("Primary generator control");
    
    fUseMeanSamplingCmd = new G4UIcmdWithABool("/primary/useMeanSampling", this);
    fUseMeanSamplingCmd->SetGuidance("Set energy sampling mode");
    fUseMeanSamplingCmd->SetGuidance("true  = Use mean energy sampling (faster)");
    fUseMeanSamplingCmd->SetGuidance("false = Use CDF-based energy sampling (more accurate)");
    fUseMeanSamplingCmd->SetParameterName("useMeanSampling", false);
    fUseMeanSamplingCmd->SetDefaultValue(true);
    fUseMeanSamplingCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    
    fAngularFileCmd = new G4UIcmdWithAString("/primary/setAngularFile", this);
    fAngularFileCmd->SetGuidance("Set the angular distribution data file path");
    fAngularFileCmd->SetParameterName("angularFile", false);
    fAngularFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    
    fMeanEnergyFileCmd = new G4UIcmdWithAString("/primary/setMeanEnergyFile", this);
    fMeanEnergyFileCmd->SetGuidance("Set the mean energy data file path");
    fMeanEnergyFileCmd->SetParameterName("meanEnergyFile", false);
    fMeanEnergyFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    
    fCDFEnergyFileCmd = new G4UIcmdWithAString("/primary/setCDFEnergyFile", this);
    fCDFEnergyFileCmd->SetGuidance("Set the CDF energy data file path");
    fCDFEnergyFileCmd->SetParameterName("cdfEnergyFile", false);
    fCDFEnergyFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fEnergyThresholdCmd = new G4UIcmdWithADoubleAndUnit("/primary/setEnergyThreshold", this);
    fEnergyThresholdCmd->SetGuidance("Set the muon energy threshold for sampling");
    fEnergyThresholdCmd->SetParameterName("energyThreshold", false);
    fEnergyThresholdCmd->SetDefaultUnit("GeV");
    fEnergyThresholdCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
    delete fUseMeanSamplingCmd;
    delete fAngularFileCmd;
    delete fMeanEnergyFileCmd;
    delete fCDFEnergyFileCmd;
    delete fEnergyThresholdCmd;
    delete fPrimaryDir;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == fUseMeanSamplingCmd) {
        fPrimaryAction->SetUseMeanSampling(fUseMeanSamplingCmd->GetNewBoolValue(newValue));
    }
    else if (command == fAngularFileCmd) {
        fPrimaryAction->SetAngularDataFile(newValue);
    }
    else if (command == fMeanEnergyFileCmd) {
        fPrimaryAction->SetMeanEnergyFile(newValue);
    }
    else if (command == fCDFEnergyFileCmd) {
        fPrimaryAction->SetCDFEnergyFile(newValue);
    }
    else if (command == fEnergyThresholdCmd) {
        fPrimaryAction->SetMuonEnergyThreshold(fEnergyThresholdCmd->GetNewDoubleValue(newValue));
    }
}

