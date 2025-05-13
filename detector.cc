#include "detector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name), fPhotonCount(0)
{
    
}

MySensitiveDetector::~MySensitiveDetector()
{}

void MySensitiveDetector::Initialize(G4HCofThisEvent*) {
    // Reset count at the beginning of each event
    fPhotonCount = 0;
}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, 
G4TouchableHistory *ROhist)
{
    // Sensitive detector (PMT front) tracks photons
    G4Track *track = aStep->GetTrack();
    G4String particleName = track->GetDefinition()->GetParticleName();
    G4double KE = track->GetKineticEnergy() / keV;
    G4int parentID = track->GetParentID();
    G4int particleID = track->GetDefinition()->GetPDGEncoding();
    G4double posZ = track->GetPosition().z() / mm;


    if ( particleName == "opticalphoton"){

        G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
        G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
        G4double time = track->GetGlobalTime() / ns;

        // Get detector positions of photons 
        const G4VTouchable* touchable = preStepPoint->GetTouchable();
        G4ThreeVector posDetector = touchable->GetTranslation();

       // Get detector volume
        G4LogicalVolume* fDetectorVolume = aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume();


        G4MaterialPropertiesTable* mpt = fDetectorVolume->GetMaterial()->GetMaterialPropertiesTable();
        G4double efficiency = 1.0;  // Default efficiency if not found in MPT
        G4double random = G4UniformRand();
        if (mpt) {
            efficiency = mpt->GetProperty("EFFICIENCY")->Value(KE *1000);
            //G4cout << "Efficiency for " << KE << " keV: " << efficiency << G4endl;
        }

        if (random < efficiency) {
            auto man = G4AnalysisManager::Instance();

            man->FillNtupleDColumn(4, 0, posDetector[0]);
            man->FillNtupleDColumn(4, 1, posDetector[1]);
            man->FillNtupleDColumn(4, 2, posDetector[2]);
            man->FillNtupleDColumn(4, 3, time);
            
            man->AddNtupleRow(4);
    
            fPhotonCount++;
        
        }
    
        //track->SetTrackStatus(fStopAndKill);
    }
    
    return true;    
}


