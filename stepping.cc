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
    if ((volume == fScoringVolume) && (particleName != "opticalphoton")){
        G4double edep = step->GetTotalEnergyDeposit();
        if (edep > 0) {
            //G4cout << "Edep: " << edep / keV << " keV in step" << step << G4endl;
            fEventAction->AddEdep(edep);
            fEventAction->AddGlobalTime(globalTime);
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

    // Store gamma scattering
    if (((particleName == "gamma"))  && (volume == fScoringVolume))
    {
        G4String processName = "None";

        const G4VProcess* postProc = step->GetPostStepPoint()->GetProcessDefinedStep();
        if (postProc) {
            processName = postProc->GetProcessName();
            //G4cout << "proc nom " << processName << G4endl;
        }
    
        // Check if it's a scattering process of interest
        if (processName == "compt" || processName == "Rayleigh" || processName == "phot") {
        // Check if this is the first interaction for this track
        if (scatteringCount[trackID] == 0) {
            G4double edep_scat = step->GetTotalEnergyDeposit();
            // Save to file
            man->FillNtupleDColumn(5, 0, posX);
            man->FillNtupleDColumn(5, 1, posY);
            man->FillNtupleDColumn(5, 2, posZ);
            man->FillNtupleDColumn(5, 3, globalTime);
            man->FillNtupleDColumn(5, 4, edep_scat);
            man->FillNtupleSColumn(5, 5, processName);
            

            // Mark this track as scattered once
            scatteringCount[trackID] = 1;

            man->FillNtupleIColumn(5, 6, scatteringCount[trackID]);
            man->AddNtupleRow(5);
        } else {
            // More than one interaction: ignore or delete if needed
            scatteringCount[trackID]++;
            return;
        }
    }

    }
    // check yields
    auto material = step->GetPreStepPoint()->GetMaterial();
    auto mpt = material->GetMaterialPropertiesTable();
    auto yield = mpt->GetProperty("SCINTILLATIONYIELD");
    if (yield) {
        G4double eDep = step->GetTotalEnergyDeposit();
        if (eDep > 0){
            G4double yieldVal = yield->Value(eDep);
            G4cout << "op photon" << G4endl;
            G4cout << "Yield at " << eDep / keV << " keV: " << yieldVal << " photons/MeV" << G4endl;
        }
    }


    // Store photon hit times at PMTs
    if ( (particleName == "opticalphoton") && (posZ > 290 || posZ < -290)) 
    {
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

        //G4cout << "parent= " << parentID << " particle id " << particleID << G4endl;
           

        // Get the creation time of the photon
        if (parentID != particleID){
    
            man->FillNtupleDColumn(3, 6, globalTime);
        }
        
        G4double edep_ph = step->GetTotalEnergyDeposit();
        man->FillNtupleDColumn(3, 7, edep_ph);

        man->AddNtupleRow(3);  

     }

     // Store photon hit times at PMTs
    if ( (particleName == "e-") && (volume == fScoringVolume)) 
    {
        // Time and true position
        man->FillNtupleDColumn(6, 0, posX);
        man->FillNtupleDColumn(6, 1, posY);
        man->FillNtupleDColumn(6, 2, posZ);
        man->FillNtupleDColumn(6, 3, globalTime);

       // G4cout << "electron = " << track->GetPosition() / mm << " mm, Time: " << globalTime << G4endl;

        // check yields
        // auto material = step->GetPreStepPoint()->GetMaterial();
        // auto mpt = material->GetMaterialPropertiesTable();
        // auto yield = mpt->GetProperty("ELECTRONSCINTILLATIONYIELD");
        // if (yield) {
        //     G4double eDep = step->GetTotalEnergyDeposit();
        //     G4double yieldVal = yield->Value(eDep);
        //    G4cout << "Yield at " << eDep / keV << " keV: " << yieldVal << " photons/MeV" << G4endl;
        // }

        
        G4double edep_e = step->GetTotalEnergyDeposit();
        man->FillNtupleDColumn(6, 4, edep_e);

        man->AddNtupleRow(6);  

     }

}
