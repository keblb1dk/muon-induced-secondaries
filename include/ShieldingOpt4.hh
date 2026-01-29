#ifndef ShieldingOpt4_h
#define ShieldingOpt4_h 1

#include "Shielding.hh"
#include "G4EmStandardPhysics_option4.hh"

class ShieldingOpt4 : public Shielding {
public:
    ShieldingOpt4() : Shielding() {
        // Replace default EM with option 4
        ReplacePhysics(new G4EmStandardPhysics_option4());
    }
};

#endif
