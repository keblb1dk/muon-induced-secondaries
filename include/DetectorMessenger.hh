#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

class DetectorMessenger : public G4UImessenger {
public:
    DetectorMessenger(DetectorConstruction*);
    virtual ~DetectorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

private:
    DetectorConstruction* fDetectorConstruction;
    G4UIdirectory* fDirectory;
    G4UIcmdWithAString* fFileNameCmd;
    G4UIcmdWithADoubleAndUnit* fRockDensityCmd;
    G4UIcmdWithAString* fRockCompositionCmd;
    G4UIcmdWithAString* fLabNameCmd;
};

#endif
