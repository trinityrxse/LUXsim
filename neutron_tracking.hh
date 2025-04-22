#ifndef NEUTRONTRACKINFO_HH
#define NEUTRONTRACKINFO_HH

#include "G4VUserTrackInformation.hh"

class NeutronTrackInfo : public G4VUserTrackInformation {
public:
    NeutronTrackInfo() : hasScattered(false) {}  // Constructor initializes flag to false
    virtual ~NeutronTrackInfo() {}

    void SetScattered() { hasScattered = true; }
    G4bool HasScattered() const { return hasScattered; }

private:
    G4bool hasScattered;
};

#endif
