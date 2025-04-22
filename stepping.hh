#ifndef STEPPING_HH
#define STEPPING_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "construction.hh"
#include "event.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "globals.hh"






class MySteppingAction : public G4UserSteppingAction
{
public:
    MySteppingAction(MyEventAction* eventAction);
    ~MySteppingAction();
    
    virtual void UserSteppingAction(const G4Step*);

private:
MyEventAction *fEventAction;

};

#endif
