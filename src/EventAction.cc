// EventAction.cc
#include "EventAction.hh"
#include "G4RunManager.hh"

EventAction::EventAction()
    : G4UserEventAction(),
    fRunAction(nullptr)
{}

EventAction::~EventAction()
{}

void EventAction::EndOfEventAction(const G4Event*)
{
    if (fRunAction) {
        const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
        fRunAction->UpdateProgress(run);
    }
}