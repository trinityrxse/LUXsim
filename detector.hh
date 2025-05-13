#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4PhysicsOrderedFreeVector.hh"  
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "construction.hh"

class MySensitiveDetector : public G4VSensitiveDetector
{
public:
    MySensitiveDetector(G4String);
    ~MySensitiveDetector();
    G4int GetPhotonCount() const { return fPhotonCount; }

private:
    void Initialize(G4HCofThisEvent*);
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *); 
    G4int fPhotonCount;

};

#endif
