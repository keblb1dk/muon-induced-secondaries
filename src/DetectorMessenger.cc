#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
    : G4UImessenger(),
    fDetectorConstruction(Det)
{

    fRockDensityCmd = new G4UIcmdWithADoubleAndUnit("/detector/setRockDensity", this);
    fRockDensityCmd->SetGuidance("Set the density of custom rock material");
    fRockDensityCmd->SetParameterName("density", false);
    fRockDensityCmd->SetUnitCategory("Volumic Mass");
    fRockDensityCmd->SetDefaultUnit("g/cm3");
    fRockDensityCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fRockCompositionCmd = new G4UIcmdWithAString("/detector/setRockComp", this);
    fRockCompositionCmd->SetGuidance("Set the rock composition from a text file");
    fRockCompositionCmd->SetGuidance("File format: Each line should contain 'ElementSymbol MassFraction'");
    fRockCompositionCmd->SetGuidance("Example: Si 0.362176");
    fRockCompositionCmd->SetParameterName("filename", false);
    fRockCompositionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fLabNameCmd = new G4UIcmdWithAString("/detector/setLabName", this);
    fLabNameCmd->SetGuidance("Set the lab name for output file naming");
    fLabNameCmd->SetGuidance("Example: snolab, lngs, soudan, hideout");
    fLabNameCmd->SetParameterName("labname", false);
    fLabNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
    delete fFileNameCmd;
    delete fRockDensityCmd;
    delete fRockCompositionCmd;
    delete fLabNameCmd;
    delete fDirectory;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == fFileNameCmd) {
        fDetectorConstruction->SetFilePattern(newValue);
    }
    else if (command == fRockDensityCmd) {
        fDetectorConstruction->SetRockDensity(fRockDensityCmd->GetNewDoubleValue(newValue));
    }
    else if (command == fRockCompositionCmd) {
        fDetectorConstruction->SetRockComposition(newValue);
    }
    else if (command == fLabNameCmd) {
        fDetectorConstruction->SetLabName(newValue);
    }
}