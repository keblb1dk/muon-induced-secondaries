#include "DetectorConstruction.hh"
#include "SteppingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4UImanager.hh"
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
#include <filesystem>
namespace fs = std::filesystem;

SteppingAction::SteppingAction(const G4String& outputDir)
    : G4UserSteppingAction(),
    fCurrentTrackID(-1),
    fOutputDir(outputDir) {
}

SteppingAction::~SteppingAction() {
}

void SteppingAction::BeginOfRunAction() {

    // Resolve lab name for output file naming
    DetectorConstruction* detector = DetectorConstruction::GetInstance();
    G4String labName = "default";
    if (detector) {
        labName = detector->GetLabName();
    }

    G4String filename = fOutputDir + "/escaped_secondary_" + labName + ".root";

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile(filename);

    G4cout << "NTuple output file opened: " << filename << G4endl;

    // Book() must be called AFTER OpenFile() so the ntuple is registered in the open file
    Book();
}

void SteppingAction::EndOfRunAction() {
    Save();
}

void SteppingAction::Save() {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
    analysisManager->Clear();
    G4cout << "\n----> NTuples are saved\n" << G4endl;
}

void SteppingAction::Initialize() {
    // This line is necessary!
    G4UImanager::GetUIpointer()->ExecuteMacroFile("simAnalysis.mac");
}

void SteppingAction::Book() {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetActivation(true);
    analysisManager->SetFirstNtupleId(0);

    // Create ntuple
    analysisManager->CreateNtuple("MuSec", "Check Outputs");

    // Column 0
    analysisManager->CreateNtupleIColumn("TrackID");
    // Column 1
    analysisManager->CreateNtupleIColumn("EventID");
    // Column 2
    analysisManager->CreateNtupleIColumn("MuonTrackID");
    // Column 3
    analysisManager->CreateNtupleDColumn("CreationEnergy (MeV)");
    // Column 4
    analysisManager->CreateNtupleDColumn("ExitEnergy (MeV)");
    // Column 5
    analysisManager->CreateNtupleDColumn("ExitX (mm)");
    // Column 6
    analysisManager->CreateNtupleDColumn("ExitY (mm)");
    // Column 7
    analysisManager->CreateNtupleDColumn("ExitZ (mm)");
    // Column 8
    analysisManager->CreateNtupleDColumn("CreationX (mm)");
    // Column 9
    analysisManager->CreateNtupleDColumn("CreationY (mm)");
    // Column 10
    analysisManager->CreateNtupleDColumn("CreationZ (mm)");
    // Column 11
    analysisManager->CreateNtupleDColumn("ExitPx");
    // Column 12
    analysisManager->CreateNtupleDColumn("ExitPy");
    // Column 13
    analysisManager->CreateNtupleDColumn("ExitPz");
    // Column 14
    analysisManager->CreateNtupleDColumn("muonPx");
    // Column 15
    analysisManager->CreateNtupleDColumn("muonPy");
    // Column 16
    analysisManager->CreateNtupleDColumn("muonPz");
    // Column 17
    analysisManager->CreateNtupleDColumn("Depth (m)");
    // Column 18
    analysisManager->CreateNtupleDColumn("TrackLength (mm)");
    // Column 19
    analysisManager->CreateNtupleDColumn("MuonTheta (deg)");
    // Column 20
    analysisManager->CreateNtupleDColumn("MuonPhi (deg)");
    // Column 21
    analysisManager->CreateNtupleDColumn("MuonInitialEnergy (GeV)");
    // Column 22
    analysisManager->CreateNtupleSColumn("creationProcess");
    // Column 23
    analysisManager->CreateNtupleSColumn("SecondaryName");

    analysisManager->FinishNtuple(0);
    analysisManager->SetNtupleActivation(true);
}


void SteppingAction::FillNTuple(const EscapedSecondaryData& data) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillNtupleIColumn(0, 0, data.trackID);
    analysisManager->FillNtupleIColumn(0, 1, data.eventID);
    analysisManager->FillNtupleIColumn(0, 2, data.muonTrackID);
    analysisManager->FillNtupleDColumn(0, 3, data.creationEnergy / MeV);
    analysisManager->FillNtupleDColumn(0, 4, data.exitEnergy / MeV);
    analysisManager->FillNtupleDColumn(0, 5, data.exitPos.x() / mm);
    analysisManager->FillNtupleDColumn(0, 6, data.exitPos.y() / mm);
    analysisManager->FillNtupleDColumn(0, 7, data.exitPos.z() / mm);
    analysisManager->FillNtupleDColumn(0, 8, data.initialPos.x() / mm);
    analysisManager->FillNtupleDColumn(0, 9, data.initialPos.y() / mm);
    analysisManager->FillNtupleDColumn(0, 10, data.initialPos.z() / mm);
    analysisManager->FillNtupleDColumn(0, 11, data.exitMomentum.x());
    analysisManager->FillNtupleDColumn(0, 12, data.exitMomentum.y());
    analysisManager->FillNtupleDColumn(0, 13, data.exitMomentum.z());
    analysisManager->FillNtupleDColumn(0, 14, data.muonMomentum.x());
    analysisManager->FillNtupleDColumn(0, 15, data.muonMomentum.y());
    analysisManager->FillNtupleDColumn(0, 16, data.muonMomentum.z());
    analysisManager->FillNtupleDColumn(0, 17, data.depth);
    analysisManager->FillNtupleDColumn(0, 18, data.trackLength / mm);
    analysisManager->FillNtupleDColumn(0, 19, data.muonTheta);
    analysisManager->FillNtupleDColumn(0, 20, data.muonPhi);
    analysisManager->FillNtupleDColumn(0, 21, data.muonInitialEnergy / GeV);
    analysisManager->FillNtupleSColumn(0, 22, data.creatorProcess);
    analysisManager->FillNtupleSColumn(0, 23, data.particleName);

    analysisManager->AddNtupleRow(0);
}


void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    if (!track) return;

    G4int trackID = track->GetTrackID();
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    // SECTION 1: Primary muon tracking
    if (track->GetParentID() == 0) {
        const G4ParticleDefinition* particle = track->GetParticleDefinition();
        if (particle == G4MuonMinus::Definition()) {
            if (track->GetCurrentStepNumber() == 1) {
                G4ThreeVector momentum = track->GetVertexMomentumDirection();
                G4ThreeVector muonMom = momentum;
                G4double cosTheta = std::abs(momentum.z());
                G4double phi = std::atan2(momentum.y(), momentum.x());
                if (phi < 0) phi += 2 * CLHEP::pi;

                MuonAngles angles;
                angles.theta = std::acos(cosTheta) * 180.0 / CLHEP::pi;
                angles.phi = phi * 180.0 / CLHEP::pi;
                angles.initialEnergy = track->GetKineticEnergy();
                angles.initialMomentum = muonMom;
                fMuonAngles[trackID] = angles;
            }
        }
        return;
    }

    // SECTION 2: Secondary tracking
    if (track->GetParentID() == 0) return;

    // Get volumes with safety checks
    G4VPhysicalVolume* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
    G4VPhysicalVolume* postVol = step->GetPostStepPoint()->GetPhysicalVolume();
    if (!preVol || !postVol) return;

    G4LogicalVolume* preVolume = preVol->GetLogicalVolume();
    G4LogicalVolume* postVolume = postVol->GetLogicalVolume();
    if (!preVolume || !postVolume) return;

    G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    // Record new secondaries created in the rock volume
    if (track->GetCurrentStepNumber() == 1 && volume->GetName() == "rockLV") {
        G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();

        SecondaryTrackInfo info;
        info.initialPos = pos;
        info.creationEnergy = track->GetVertexKineticEnergy();

        G4int parentID = track->GetParentID();
        if (fMuonAngles.find(parentID) != fMuonAngles.end()) {
            // Direct child of a primary muon
            info.muonTrackID = parentID;
        }
        else {
            // Higher-generation secondary — inherit ancestor muon from parent's record
            auto parentIter = fSecondaryTracks.find(parentID);
            info.muonTrackID = (parentIter != fSecondaryTracks.end())
                ? parentIter->second.muonTrackID
                : -1;
        }

        fSecondaryTracks[trackID] = info;

        const G4VProcess* creatorProcess = track->GetCreatorProcess();
    }

    // Process existing secondary tracks
    auto trackIter = fSecondaryTracks.find(trackID);
    if (trackIter == fSecondaryTracks.end()) return;

    G4LogicalVolume* nextVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    // Check for volume boundary crossing into lab
    if (nextVolume && step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary && nextVolume->GetName() == "labLV") {
        G4ThreeVector crossingPoint = step->GetPreStepPoint()->GetPosition();
        G4double exitEnergy = step->GetPreStepPoint()->GetKineticEnergy();
        G4ThreeVector exitMomentum = step->GetPreStepPoint()->GetMomentumDirection();

        const G4VProcess* creatorProcess = track->GetCreatorProcess();

        // Define lab dimensions
        G4double labX = 3.05 / 2.0;
        G4double labY = 5.2 / 2.0;
        G4double labZ = 3.05;

        try {
            const G4ThreeVector& initialPos = trackIter->second.initialPos;
            G4double trackLength = std::sqrt(
                std::pow(crossingPoint.x() - initialPos.x(), 2) +
                std::pow(crossingPoint.y() - initialPos.y(), 2) +
                std::pow(crossingPoint.z() - initialPos.z(), 2));

            G4double depth = 0.0;

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

            // Get muon info
            G4double muonTheta = -1;
            G4double muonPhi = -1;
            G4double muonInitialEnergy = -1;
            G4ThreeVector muonMomentum;
            auto muonIter = fMuonAngles.find(trackIter->second.muonTrackID);
            if (muonIter != fMuonAngles.end()) {
                muonTheta = muonIter->second.theta;
                muonPhi = muonIter->second.phi;
                muonInitialEnergy = muonIter->second.initialEnergy;
                muonMomentum = muonIter->second.initialMomentum;
            }

            // Build data struct and fill NTuple directly
            EscapedSecondaryData data;
            data.trackID = trackID;
            data.eventID = eventID;
            data.muonTrackID = trackIter->second.muonTrackID;
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

            FillNTuple(data);

            fSecondaryTracks.erase(trackID);
        }
        catch (const std::exception& e) {
            G4cerr << "Error processing secondary escape: " << e.what() << G4endl;
        }

        track->SetTrackStatus(fStopAndKill);
    }
}

void SteppingAction::Reset() {
    fSecondaryTracks.clear();
    fMuonAngles.clear();
}

void SteppingAction::EndOfEventAction() {
    Reset();
}
