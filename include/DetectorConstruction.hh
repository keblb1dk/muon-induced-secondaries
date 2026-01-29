#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include <map>
#include <string>

class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
	
	    // Methods to set custom rock properties
    void SetRockDensity(G4double density);
    void SetRockComposition(const G4String& filename);

    // Lab name for output files
    void SetLabName(const G4String& name);
    G4String GetLabName() const { return fLabName; }
    static DetectorConstruction* GetInstance() { return fInstance; }

    void SetFilePattern(const G4String& pattern) { fFilePattern = pattern; }
	
private:
    void LoadRockComposition(const G4String& filename);
    void PrintRockComposition();
    G4Material* CreateCustomRock();

    DetectorMessenger* fMessenger;

    // Static instance for global access
    static DetectorConstruction* fInstance;

    // Custom rock properties
    G4double fCustomRockDensity;
    std::map<std::string, G4double> fRockComposition;  // element symbol -> mass fraction
    G4bool fUseCustomRock;

    // Lab name for output files
    G4String fLabName;

    // Existing properties
    G4String fFilePattern;
};

#endif
