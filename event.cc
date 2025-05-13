#include "event.hh"
#include <cfloat> // For DBL_MAX

MyEventAction::MyEventAction(MyRunAction*)
    : fEdep(0.), fTime(DBL_MAX), fPhotons(0), fPhotonsPhotocath(0)
{
    // Constructor body left empty intentionally
}

MyEventAction::~MyEventAction() {}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
    fEdep = 0.;
    fTime = DBL_MAX;
    fPhotons = 0;
    fPhotonsPhotocath = 0;
    fPhotonTrackIDs.clear();
    fPhotonDetTrackIDs.clear();
}

void MyEventAction::EndOfEventAction(const G4Event*)
{   
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    G4VSensitiveDetector* sd = sdManager->FindSensitiveDetector("SensitiveDetector");

    if (!sd) {
        G4cerr << "Sensitive detector not found!" << G4endl;
        return;
    }

    MySensitiveDetector* mySD = dynamic_cast<MySensitiveDetector*>(sd);
    if (!mySD) {
        G4cerr << "Failed to cast to MySensitiveDetector!" << G4endl;
        return;
    }

    G4int count = mySD->GetPhotonCount();
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    // Total energy deposited in xenon over the whole event
    man->FillNtupleDColumn(1, 0, fEdep);
    man->FillNtupleDColumn(1, 1, fTime);
    man->FillNtupleDColumn(1, 2, fPhotons);
    man->FillNtupleDColumn(1, 3, fPhotonsPhotocath);
    man->FillNtupleDColumn(1, 4, count); // Photons Detected
    man->AddNtupleRow(1);
}
