#include "SteppingAction.hh"
#include "RunAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4Neutron.hh"
#include "G4MuonMinus.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include <vector>
#include <tuple>
#include <unordered_map>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;
#include <iomanip>

SteppingAction::SteppingAction(const G4String& outputDir)
    : G4UserSteppingAction(),
    fCurrentTrackID(-1),
    fOutputDir(outputDir),
    bufferSize(10000) {  // Adjust buffer size as needed

    std::filesystem::path dirPath(fOutputDir.data());
    std::filesystem::create_directories(dirPath);

    size_t pos = fOutputDir.find_last_of('_');
    G4String coreNum = "";
    if (pos != G4String::npos) {
        coreNum = "_" + fOutputDir.substr(pos + 1);
    }

    // Construct the file path
    G4String filePath = fOutputDir + "/escaped_secondaries.csv";
    fEscapedSecondaryData.open(filePath, std::ios::out | std::ios::app);

    if (!fEscapedSecondaryData.is_open()) {
        G4Exception("SteppingAction::SteppingAction", "OutputError",
            FatalException,
            G4String("Could not open output file: " + filePath).c_str());
    }
    // write header of csv if empty
    if (fEscapedSecondaryData.tellp() == 0) {
        fEscapedSecondaryData << "TrackID,EventID,CreationEnergy (MeV),ExitEnergy (MeV),"
            << "ExitX (mm),ExitY (mm),ExitZ (mm),"
            << "CreationX (mm),CreationY (mm),CreationZ (mm),"
            << "ExitPx,ExitPy,ExitPz,"
            << "muonPx,muonPy,muonPz,"
            << "Depth (m),"
            << "TrackLength (mm),MuonTheta (deg),MuonPhi (deg),MuonInitialEnergy (GeV),creationProcess,secondaryName" << std::endl;
    }
}

SteppingAction::~SteppingAction() {
    FlushBuffer();  // Ensure any remaining data is written
    if (fEscapedSecondaryData.is_open()) {
        fEscapedSecondaryData.close();
    }
}

void SteppingAction::FlushBuffer() {
    // write data to csv
    if (fEscapedSecondaryData.is_open()) {
        for (const auto& data : dataBuffer) {
            fEscapedSecondaryData << std::fixed << std::setprecision(15)
                << data.trackID << ","
                << data.eventID << ","
                << data.creationEnergy / MeV << ","
                << data.exitEnergy / MeV << ","
                << data.exitPos.x() / mm << ","
                << data.exitPos.y() / mm << ","
                << data.exitPos.z() / mm << ","
                << data.initialPos.x() / mm << ","
                << data.initialPos.y() / mm << ","
                << data.initialPos.z() / mm << ","
                << data.exitMomentum.x() << ","
                << data.exitMomentum.y() << ","
                << data.exitMomentum.z() << ","
                << data.muonMomentum.x() << ","
                << data.muonMomentum.y() << ","
                << data.muonMomentum.z() << ","
                << data.depth << ","
                << data.trackLength / mm << ","
                << data.muonTheta << ","
                << data.muonPhi << ","
                << data.muonInitialEnergy / GeV << ","
                << data.creatorProcess << ","
                << data.particleName
                << std::endl;
        }
    }
    dataBuffer.clear();
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    if (!track) return;

    G4int trackID = track->GetTrackID();
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    // SECTION 1: Primary muon tracking
    if (trackID == 1) {
        const G4ParticleDefinition* particle = track->GetParticleDefinition();
        if (particle == G4MuonMinus::Definition()) {
            if (track->GetCurrentStepNumber() == 1) {

                // get primary muon information
                G4ThreeVector momentum = track->GetVertexMomentumDirection();
                G4ThreeVector muonMom = momentum;
                G4double cosTheta = std::abs(momentum.z());
                G4double phi = std::atan2(momentum.y(), momentum.x());
                if (phi < 0) phi += 2 * CLHEP::pi; // ensure phi stays 0 - 360

                MuonAngles angles;
                angles.theta = std::acos(cosTheta) * 180.0 / CLHEP::pi;
                angles.phi = phi * 180.0 / CLHEP::pi;
                angles.initialEnergy = track->GetKineticEnergy();
                angles.initialMomentum = muonMom;
                fMuonAngles[eventID] = angles;
            }
        }
        return;
    }

    const G4ParticleDefinition* particle = track->GetParticleDefinition();

    // SECTION 2: Secondary tracking
    if (track->GetTrackID() == 1) return;  // ignore primary muons
    

    // Get volumes with safety checks
    G4VPhysicalVolume* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
    G4VPhysicalVolume* postVol = step->GetPostStepPoint()->GetPhysicalVolume();
    if (!preVol || !postVol) return;

    G4LogicalVolume* preVolume = preVol->GetLogicalVolume();
    G4LogicalVolume* postVolume = postVol->GetLogicalVolume();
    if (!preVolume || !postVolume) return;

    G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    // Record new secondaries with boundary conditions
    // must be created in the rock volume
    if (track->GetCurrentStepNumber() == 1 && volume->GetName() == "rShell") {
        G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();

        SecondaryTrackInfo info;
        info.initialPos = pos;
        info.creationEnergy = track->GetVertexKineticEnergy();
        fSecondaryTracks[trackID] = info;

        const G4VProcess* creatorProcess = track->GetCreatorProcess();
    }

    // Process existing secondary tracks
    auto trackIter = fSecondaryTracks.find(trackID);
    if (trackIter == fSecondaryTracks.end()) return; // Only process tracked secondaries

    G4LogicalVolume* nextVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    // Check for volume boundary crossing into lab (only checks first entry and ignores backscattering)
    if (nextVolume && step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary && nextVolume->GetName() == "labLV") {
        G4ThreeVector prePos = step->GetPreStepPoint()->GetPosition();
        G4ThreeVector postPos = step->GetPostStepPoint()->GetPosition();

        G4ThreeVector crossingPoint = prePos;
        G4double exitEnergy = step->GetPreStepPoint()->GetKineticEnergy();
        G4ThreeVector exitMomentum = step->GetPreStepPoint()->GetMomentumDirection();

        const G4VProcess* creatorProcess = track->GetCreatorProcess();

        //define lab dimensions
        G4double labX = 3.0 / 2.0;
        G4double labY = 5.0 / 2.0;
        G4double labZ = 3.0;

        try {
            const G4ThreeVector& initialPos = trackIter->second.initialPos;
            G4double trackLength = std::sqrt(
                std::pow(crossingPoint.x() - initialPos.x(), 2) +
                std::pow(crossingPoint.y() - initialPos.y(), 2) +
                std::pow(crossingPoint.z() - initialPos.z(), 2));

            G4double depth = 0.0;
            // straight-line distance from lab calculations
                // for z <= labZ …. Below the top of the lab
                // 1D distance from x faces
            if ((labX * m <= initialPos.x() || initialPos.x() <= -labX * m) &&
                (-labY * m <= initialPos.y() && initialPos.y() <= labY * m) &&
                (initialPos.z() <= labZ * m)) {
                depth = abs(initialPos.x() / m) - labX;
            }
            // 1D distance from y faces
            else if ((labY * m <= initialPos.y() || initialPos.y() <= -labY * m) &&
                (-labX * m <= initialPos.x() && initialPos.x() <= labX * m) &&
                (initialPos.z() <= labZ * m)) {
                depth = abs(initialPos.y() / m) - labY;
            }
            // 1D distance from z face
            else if ((-labY * m <= initialPos.y() && initialPos.y() <= labY * m) &&
                (-labX * m <= initialPos.x() && initialPos.x() <= labX * m) &&
                (labZ * m <= initialPos.z())) {
                depth = initialPos.z() / m - labZ;
            }
            // 2D distance from x > labX/2 and y > labY/2
            else if ((labX * m < initialPos.x()) && (labY * m < initialPos.y()) &&
                (initialPos.z() <= labZ * m)) {
                depth = std::sqrt(pow(abs(initialPos.x()) / m - labX, 2) + pow(abs(initialPos.y()) / m - labY, 2));
            }
            // 2D distance from x > labX/2 and y < -labY/2
            else if ((labX * m < initialPos.x()) && (-labY * m > initialPos.y()) &&
                (initialPos.z() <= labZ * m)) {
                depth = std::sqrt(pow(abs(initialPos.x()) / m - labX, 2) + pow(abs(initialPos.y()) / m - labY, 2));
            }
            // 2D distance from x < -labX/2 and y < -labY/2
            else if ((-labX * m > initialPos.x()) && (-labY * m > initialPos.y()) &&
                (initialPos.z() <= labZ * m)) {
                depth = std::sqrt(pow(abs(initialPos.x()) / m - labX, 2) + pow(abs(initialPos.y()) / m - labY, 2));
            }
            // 2D distance from x < -labX/2 and y > labY/2
            else if ((-labX * m > initialPos.x()) && (labY * m < initialPos.y()) &&
                (initialPos.z() <= labZ * m)) {
                depth = std::sqrt(pow(abs(initialPos.x()) / m - labX, 2) + pow(abs(initialPos.y()) / m - labY, 2));
            }

            // for z > labZ ... above the top of the lab
            // 2D distance from top of x faces
            else if ((labX * m <= initialPos.x() || initialPos.x() <= -labX * m) &&
                (-labY * m <= initialPos.y() && initialPos.y() <= labY * m) &&
                (initialPos.z() > labZ * m)) {
                depth = std::sqrt(pow(abs(initialPos.x()) / m - labX, 2) + pow(initialPos.z() / m - labZ, 2));
            }
            // 2D distance from top of y faces
            else if ((labY * m <= initialPos.y() || initialPos.y() <= -labY * m) &&
                (-labX * m <= initialPos.x() && initialPos.x() <= labX * m) &&
                (initialPos.z() > labZ * m)) {
                depth = std::sqrt(pow(abs(initialPos.y()) / m - labY, 2) + pow(initialPos.z() / m - labZ, 2));
            }
            // 3D distance from x > 1.5 and y > 1.5
            else if ((labX * m < initialPos.x()) && (labY * m < initialPos.y()) &&
                (initialPos.z() > labZ * m)) {
                depth = std::sqrt(pow(abs(initialPos.x()) / m - labX, 2) +
                    pow(abs(initialPos.y()) / m - labY, 2) +
                    pow(initialPos.z() / m - labZ, 2));
            }
            // 3D distance from x > 1.5 and y < -1.5
            else if ((labX * m < initialPos.x()) && (-labY * m > initialPos.y()) &&
                (initialPos.z() > labZ * m)) {
                depth = std::sqrt(pow(abs(initialPos.x()) / m - labX, 2) +
                    pow(abs(initialPos.y()) / m - labY, 2) +
                    pow(initialPos.z() / m - labZ, 2));
            }
            // 3D distance from x < -1.5 and y < -1.5
            else if ((-labX * m > initialPos.x()) && (-labY * m > initialPos.y()) &&
                (initialPos.z() > labZ * m)) {
                depth = std::sqrt(pow(abs(initialPos.x()) / m - labX, 2) +
                    pow(abs(initialPos.y()) / m - labY, 2) +
                    pow(initialPos.z() / m - labZ, 2));
            }
            // 3D distance from x < -1.5 and y > 1.5
            else if ((-labX * m > initialPos.x()) && (labY * m < initialPos.y()) &&
                (initialPos.z() > labZ * m)) {
                depth = std::sqrt(pow(abs(initialPos.x()) / m - labX, 2) +
                    pow(abs(initialPos.y()) / m - labY, 2) +
                    pow(initialPos.z() / m - labZ, 2));
            }

            else {
                G4cerr << "Warning: No valid depth calculation case for position ("
                    << initialPos.x() / m << ", " << initialPos.y() / m << ", "
                    << initialPos.z() / m << ")" << G4endl;
            }

            // Get muon info after it created a secondary
            G4double muonTheta = -1;
            G4double muonPhi = -1;
            G4double muonInitialEnergy = -1;
            G4ThreeVector muonMomentum;
            auto muonIter = fMuonAngles.find(eventID);
            if (muonIter != fMuonAngles.end()) {
                muonTheta = muonIter->second.theta;
                muonPhi = muonIter->second.phi;
                muonInitialEnergy = muonIter->second.initialEnergy;
                muonMomentum = muonIter->second.initialMomentum;
            }

            // Add data to buffer
            EscapedSecondaryData data;
            data.trackID = trackID;
            data.eventID = eventID;
            data.creationEnergy = track->GetVertexKineticEnergy();
            data.exitEnergy = exitEnergy;
            data.exitPos = crossingPoint;
            data.initialPos = initialPos;
            data.exitMomentum = exitMomentum;
            data.depth = depth;
            data.trackLength = trackLength;
            data.muonTheta = muonTheta;
            data.muonPhi = muonPhi;
            data.muonInitialEnergy = muonInitialEnergy;
            data.muonMomentum = muonMomentum;
            data.creatorProcess = creatorProcess ? creatorProcess->GetProcessName() : "Unknown";
            data.particleName = track->GetParticleDefinition()->GetParticleName();

            dataBuffer.push_back(data);

            // Flush buffer if it reaches the specified size
            if (dataBuffer.size() >= bufferSize) {
                FlushBuffer();
            }

            // Update count and cleanup
            fSecondaryTracks.erase(trackID);  // Remove after successful processing
        }
        catch (const std::exception& e) {
            G4cerr << "Error processing secondary escape: " << e.what() << G4endl;
        }

        track->SetTrackStatus(fStopAndKill); // kill the secondary track when we're done with it
    }
}
