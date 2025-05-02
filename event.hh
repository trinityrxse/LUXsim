#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "run.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryVertex.hh"
#include "globals.hh"

class MyEventAction : public G4UserEventAction
{
public: 
MyEventAction(MyRunAction*);
~MyEventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void AddEdep(G4double edep) {  fEdep += edep;  }
    void AddGlobalTime(G4double globalTime)
    {
        // we want to save the start time of the energy deposition 
        if (globalTime < fTime)
        {fTime = globalTime;}
    }
private:
    G4double fEdep;
    G4double fTime;
};
#endif