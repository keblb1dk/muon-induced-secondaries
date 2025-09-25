#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH
#include <fstream>

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);
private:
    G4ParticleGun* fParticleGun;
    std::ofstream fTotalAngles;
};

#endif
