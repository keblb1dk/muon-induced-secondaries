// EventAction.hh
#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "RunAction.hh"

class EventAction : public G4UserEventAction {
public:
    EventAction();  // Default constructor
    virtual ~EventAction();

    void SetRunAction(RunAction* runAction) { fRunAction = runAction; }
    virtual void EndOfEventAction(const G4Event*);

private:
    RunAction* fRunAction;
};

#endif