#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <iostream>
#include "G4ios.hh"

RunAction::RunAction()
    : G4UserRunAction(),
    fPreviousEventCount(0) {
    
}

RunAction::~RunAction() {
  
}

void RunAction::BeginOfRunAction(const G4Run*) {

    fPreviousEventCount = 0;
    fStartTime = std::chrono::high_resolution_clock::now();
    G4cout << "Starting new run..." << G4endl;
}

void RunAction::EndOfRunAction(const G4Run* run) {
    G4int nEvents = run->GetNumberOfEvent();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = endTime - fStartTime;
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);

    G4cout << "   \n========== Final Run Summary ==========\n"
        << "   Total Events Completed: " << nEvents << " / " << run->GetNumberOfEventToBeProcessed()
        << "    (100.00%)\n"
        << "   Final Time: " << minutes.count() << " minutes "
        << (seconds.count() % 60) << " seconds\n"
        << "   Total seconds: " << seconds.count() << "\n"
        << "   =====================================" << G4endl;

}

// print progress every N events
void RunAction::UpdateProgress(const G4Run* run) {
    G4int nEvents = run->GetNumberOfEvent();
    if (nEvents - fPreviousEventCount >= 20000) {
        auto currentTime = std::chrono::high_resolution_clock::now();
        auto duration = currentTime - fStartTime;
        auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
        auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);

        G4double percentage = 100.0 * nEvents / run->GetNumberOfEventToBeProcessed();

        G4cout << "   \n=== Progress Update ===\n"
            << "   Events: " << nEvents << " / " << run->GetNumberOfEventToBeProcessed()
            << "    (" << std::fixed << std::setprecision(2) << percentage << "%)\n"
            << "   Time elapsed: " << minutes.count() << " minutes "
            << (seconds.count() % 60) << " seconds\n"
            << "   ======================" << G4endl;

        fPreviousEventCount = nEvents;
    }
}