#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include <unordered_map>
#include <fstream>
#include <vector>
#include <string>

class G4Step;

// Define a structure to hold the data temporarily
struct EscapedSecondaryData {
    G4int trackID;
    G4int eventID;
    G4double creationEnergy;
    G4double exitEnergy;
    G4ThreeVector exitPos;
    G4ThreeVector initialPos;
    G4ThreeVector exitMomentum;
    G4ThreeVector muonMomentum;
    G4double depth;
    G4double trackLength;
    G4double muonTheta;
    G4double muonPhi;
    G4double muonInitialEnergy;
    std::string creatorProcess;
    std::string particleName;
};

struct MuonAngles {
    G4double theta;
    G4double phi;
    G4double initialEnergy;
    G4ThreeVector initialMomentum;
};

struct SecondaryTrackInfo {
    G4ThreeVector initialPos;
    G4double creationEnergy = 0.0;
};

class SteppingAction : public G4UserSteppingAction {
public:

    SteppingAction(const G4String& outputDir);
    virtual ~SteppingAction();
    virtual void UserSteppingAction(const G4Step*);

private:
    G4String fOutputDir;

    G4int fCurrentTrackID;
    // create maps for the data
    std::unordered_map<G4int, MuonAngles> fMuonAngles;
    std::unordered_map<G4int, SecondaryTrackInfo> fSecondaryTracks;
    std::ofstream fEscapedSecondaryData;
    G4double depth = 0.0;

    // Buffered data storage
    std::vector<EscapedSecondaryData> dataBuffer;
    size_t bufferSize; // Size of the buffer

    // Method to flush the buffer to the file
    void FlushBuffer();
};

#endif