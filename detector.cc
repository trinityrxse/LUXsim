#include "detector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
    
}

MySensitiveDetector::~MySensitiveDetector()
{}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, 
G4TouchableHistory *ROhist)
{
    // Sensitive detector (PMT front) tracks photons
    G4Track *track = aStep->GetTrack();
    G4String particleName = track->GetDefinition()->GetParticleName();
    G4int parentID = track->GetParentID();
    G4int particleID = track->GetDefinition()->GetPDGEncoding();
    G4double posZ = track->GetPosition().z() / mm;

    //G4cout << "nom: " << particleName << G4endl;
    // G4cout << "Z: " << posZ << G4endl;
    //G4cout << "actual pos: " << track->GetPosition() / mm << G4endl;

    if ( (particleName == "opticalphoton") && (posZ > 274 || posZ < -274)){

        G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
        G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
        G4double time = track->GetGlobalTime() / ns;

        // Get detector positions of photons 
        const G4VTouchable* touchable = preStepPoint->GetTouchable();
        G4ThreeVector posDetector = touchable->GetTranslation();

       // G4cout << "parent= " << parentID << " particle id " << particleID << G4endl;
           

        auto man = G4AnalysisManager::Instance();

        man->FillNtupleDColumn(4, 0, posDetector[0]);
        man->FillNtupleDColumn(4, 1, posDetector[1]);
        man->FillNtupleDColumn(4, 2, posDetector[2]);
        man->FillNtupleDColumn(4, 3, time);

        // Get the creation time of the photon
        if (parentID != particleID){
            man->FillNtupleDColumn(4, 4, time);
        }
        
        man->AddNtupleRow(4);

        //track->SetTrackStatus(fStopAndKill);
    }
    
    return true;

    
}


