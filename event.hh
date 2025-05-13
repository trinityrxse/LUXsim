#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "run.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryVertex.hh"
#include "globals.hh"
#include <set>
#include "G4SDManager.hh"
#include "detector.hh"

class MyEventAction : public G4UserEventAction
{
public: 
    MyEventAction(MyRunAction*);
    ~MyEventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void AddEdep(G4double edep) { fEdep += edep; }

    void AddGlobalTime(G4double globalTime)
    {
        // Store the earliest time of energy deposition
        if (globalTime < fTime) {
            fTime = globalTime;
        }
    }

    void AddPhotonTrackID(G4int trackID) {
        if (fPhotonTrackIDs.insert(trackID).second) {
            fPhotons += 1;
        }
    }

    void AddPhotonDetTrackID(G4int trackID) {
            fPhotonsPhotocath += 1;
    }

private:
    G4double fEdep = 0;
    G4double fTime = DBL_MAX;
    G4int fPhotons = 0;
    std::set<G4int> fPhotonTrackIDs;
    G4int fPhotonsPhotocath = 0;
    std::set<G4int> fPhotonDetTrackIDs;
};

#endif
