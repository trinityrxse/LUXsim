#include "stepping.hh"
#include "neutron_tracking.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
    fEventAction = eventAction;
}
    

MySteppingAction::~MySteppingAction()
{
}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{

    G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    // Ensure scoring volume is retrieved safely
    const MyDetectorConstruction *detectorConstruction = 
        static_cast<const MyDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    // Scoring volume = xenon
    G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();
    if (!fScoringVolume) {
        G4cerr << "Error: Scoring volume not defined!" << G4endl;
        return;
    }
    

    // Get the particle track
    G4Track* track = step->GetTrack();
    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4StepPoint* postStepPoint = step->GetPostStepPoint();
    G4String particleName = track->GetDefinition()->GetParticleName();
    G4int particleID = track->GetDefinition()->GetPDGEncoding();
    G4double time = track->GetGlobalTime() / ns;  // Convert to nanoseconds
    G4String volumeName = track->GetVolume()->GetName();

    // ROOT File Storage
    auto man = G4AnalysisManager::Instance();

    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4int trackID = track->GetTrackID();
    G4int parentID = track->GetParentID();
    G4double globalTime = track->GetGlobalTime() / ns;
    G4double KE = track->GetKineticEnergy() / keV;
    G4double posX = track->GetPosition().x() / mm;
    G4double posY = track->GetPosition().y() / mm;
    G4double posZ = track->GetPosition().z() / mm;

    // All particles in event information at first and last step
    // Store particle entry/exit in Xenon
    if (volume == fScoringVolume) {
        if (step->IsFirstStepInVolume() || step->IsLastStepInVolume()) {
            man->FillNtupleIColumn(0, 0, eventID);   
            man->FillNtupleIColumn(0, 1, trackID);   
            man->FillNtupleIColumn(0, 2, parentID);  
            man->FillNtupleIColumn(0, 3, particleID); 
            man->FillNtupleDColumn(0, 4, globalTime);
            man->FillNtupleDColumn(0, 5, KE);        
            man->FillNtupleDColumn(0, 6, posX);       
            man->FillNtupleDColumn(0, 7, posY);       
            man->FillNtupleDColumn(0, 8, posZ);       
            man->AddNtupleRow(0);  
        }
    }

    // Store energy deposition in scoring volume (xenon)
    if (volume == fScoringVolume) {
        G4double edep = step->GetTotalEnergyDeposit();
        if (edep > 0) {
            //G4cout << "Edep: " << edep / keV << " keV in step" << step << G4endl;
            fEventAction->AddEdep(edep / keV);
        }
    }

    // Check if it's a neutron in the xenon
    if ((particleName == "neutron") && (volume == fScoringVolume) ){

        // Retrieve or attach NeutronTrackInfo
        NeutronTrackInfo* trackInfo = dynamic_cast<NeutronTrackInfo*>(track->GetUserInformation());
        if (!trackInfo) {
            trackInfo = new NeutronTrackInfo();
            track->SetUserInformation(trackInfo);
        }

        // If the neutron has already scattered, ignore further events
        if (trackInfo->HasScattered()) return;

        // Get the process name
        const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();

        if (process) {
            G4String processName = process->GetProcessName();

            // Check if the process is elastic scattering
            if (processName == "hadElastic") {

                // Get neutron initial and final momenta
                G4ThreeVector p_initial = step->GetPreStepPoint()->GetMomentum();
                G4ThreeVector p_final = step->GetPostStepPoint()->GetMomentum();

                const G4StepPoint* postStep = step->GetPostStepPoint();
                if (postStep) {
                    G4Material* mat = postStep->GetMaterial();
                    if (mat && mat->GetNumberOfElements() > 0) {
                        G4double nucleusMass = mat->GetElement(0)->GetAtomicMassAmu() * CLHEP::amu_c2;

                        // Compute recoil momentum
                        G4ThreeVector p_recoil = p_initial - p_final;
                        G4double p_recoil_mag = p_recoil.mag();

                        // Compute recoil energy: E_recoil = p^2 / (2M)
                        G4double E_recoil = (std::pow(p_recoil_mag, 2) / (2 * nucleusMass)) ; //in MeV

                        // Get the momentum before and after scattering
                        G4ThreeVector p_initial_dir = step->GetPreStepPoint()->GetMomentumDirection();
                        G4ThreeVector p_final_dir = step->GetPostStepPoint()->GetMomentumDirection();

                        // Compute the scattering angle
                        G4double scatteringAngle_rad = p_initial_dir.angle(p_final_dir);

                        // Convert to degrees
                        G4double scatteringAngle = scatteringAngle_rad * (180. / CLHEP::pi);

                        //G4cout << "Scattering angle: " << scatteringAngle << " degrees" << G4endl;

                        // Find energy deposited in xenon
                        G4double neutronEdep = step->GetTotalEnergyDeposit();
                        //G4cout << " neutron edep" << neutronEdep /keV << " keV" << G4endl;
                        // check that this matches (or closely resembles) recoil energy

                        // CHECK: For elastic collision, E_neutron_initial = E_neutron_final + E_recoil
                        G4double E_initial = step->GetPreStepPoint()->GetKineticEnergy();
                        G4double E_final = step->GetPostStepPoint()->GetKineticEnergy();
                        G4double E_elastic = E_final + E_recoil;
                        //G4double E_elastic_SA = E_final + E_recoil_SA;
                        //G4cout << "E_neutron_final + E_recoil: " << E_elastic / keV << " keV, Initial:" << E_initial / keV << " keV" << G4endl;
            
                        man->FillNtupleDColumn(2, 0, scatteringAngle);   
                        man->FillNtupleDColumn(2, 1, posX);   
                        man->FillNtupleDColumn(2, 2, posY);  
                        man->FillNtupleDColumn(2, 3, posZ); 
                        man->FillNtupleDColumn(2, 4, E_recoil);
                        man->FillNtupleDColumn(2, 5, neutronEdep);        
                        man->FillNtupleDColumn(2, 6, globalTime);            
                        man->AddNtupleRow(2);  
                    }
                }
                

                // Mark this neutron so it does not register another scatter
                trackInfo->SetScattered();

            }
        }
    }


    // G4cout << "posZ " << posZ << ", particle type: " << particleName << G4endl;

    // Store photon hit times at PMTs
    if ( (particleName == "opticalphoton" || particleName == "gamma") && (posZ > 274 || posZ < -274)) {

        // Time and true position
        man->FillNtupleDColumn(3, 0, globalTime);
        man->FillNtupleDColumn(3, 1, posX);
        man->FillNtupleDColumn(3, 2, posY);
        man->FillNtupleDColumn(3, 3, posZ);

        //G4cout << "Optical PhotonPosition = " << track->GetPosition() / mm << " mm, Time: " << globalTime << G4endl;

        // Store optical photon energy and wavelength if it is detected by PMTs
        G4double energy = track->GetTotalEnergy();
        if (energy > 0){
            G4double wavelength = (1240.0 / (energy / eV));  // Convert to nm
            //G4cout << "Optical Photon: Energy = " << energy / keV << " keV, Wavelength (nm): " << wavelength << G4endl;
            man->FillNtupleDColumn(3, 4, energy);
            man->FillNtupleDColumn(3, 5, wavelength);
        }

        man->AddNtupleRow(3);  

     }
}
