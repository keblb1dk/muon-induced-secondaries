#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <chrono>
#include <fstream>

class G4Run;

class RunAction : public G4UserRunAction {
public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    void UpdateProgress(const G4Run* run);


private:
    std::chrono::high_resolution_clock::time_point fStartTime;

    G4int fPreviousEventCount;
};

#endif