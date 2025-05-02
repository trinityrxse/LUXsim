#include "event.hh"

MyEventAction::MyEventAction(MyRunAction*)
{

}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
    fEdep = 0.;
    fTime = 9999999999999.;
}

void MyEventAction::EndOfEventAction(const G4Event*)
{

    G4AnalysisManager *man = G4AnalysisManager::Instance();


    // Total energy deposited in xenon over the whole event
    man->FillNtupleDColumn(1, 0, fEdep);
    man->FillNtupleDColumn(1, 1, fTime);
    man->AddNtupleRow(1);



}