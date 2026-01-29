#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <fstream>

class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

    // Methods to control energy sampling mode
    void SetUseMeanSampling(G4bool val);
    G4bool GetUseMeanSampling() const { return fUseMeanSampling; }

    //MEthod to set muon energy threshold...default 0.0
    void SetMuonEnergyThreshold(G4double eThreshold);
    G4double GetMuonEnergyThreshold() const { return fEnergyThreshold; }


    // Methods to set data file paths
    void SetAngularDataFile(const G4String& filename);
    void SetMeanEnergyFile(const G4String& filename);
    void SetCDFEnergyFile(const G4String& filename);

private:
    G4ParticleGun* fParticleGun;
    std::ofstream fTotalAngles;
    PrimaryGeneratorMessenger* fMessenger;

    // Control flag for energy sampling mode
    G4bool fUseMeanSampling;  // true = mean energy, false = CDF sampling
    G4double fEnergyThreshold;

    // File paths
    G4String fAngularDataFile;
    G4String fMeanEnergyFile;
    G4String fCDFEnergyFile;
};

#endif
