#include "run.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"

MyRunAction::MyRunAction() {
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    man->CreateNtuple("Hits", "Hits");
    man->CreateNtupleIColumn("eventID");
    man->CreateNtupleIColumn("trackID");
    man->CreateNtupleIColumn("parentID");
    man->CreateNtupleIColumn("particleID");
    man->CreateNtupleDColumn("globalTime");
    man->CreateNtupleDColumn("Kinetic Energy");
    man->CreateNtupleDColumn("posX");
    man->CreateNtupleDColumn("posY");
    man->CreateNtupleDColumn("posZ");
    man->FinishNtuple(0);

    man->CreateNtuple("Scoring", "Scoring");
    man->CreateNtupleDColumn("fEdep");
    man->CreateNtupleDColumn("startTime");
    man->FinishNtuple(1);


    man->CreateNtuple("Neutrons", "Neutrons");
    man->CreateNtupleDColumn("elasticScatteringAngles");
    man->CreateNtupleDColumn("posX");
    man->CreateNtupleDColumn("posY");
    man->CreateNtupleDColumn("posZ");
    man->CreateNtupleDColumn("Recoil Energy");
    man->CreateNtupleDColumn("Energy Deposited");
    man->CreateNtupleDColumn("Time");
    man->FinishNtuple(2);

    man->CreateNtuple("Photons", "Photons");
    man->CreateNtupleDColumn("photonTimes");
    man->CreateNtupleDColumn("posX");
    man->CreateNtupleDColumn("posY");
    man->CreateNtupleDColumn("posZ");
    man->CreateNtupleDColumn("energies");
    man->CreateNtupleDColumn("wavelengths");
    man->CreateNtupleDColumn("creationTime");
    man->CreateNtupleDColumn("edep");
    man->FinishNtuple(3);

    man->CreateNtuple("Photon Detector Hits", "Photon Detector Hits");
    man->CreateNtupleDColumn("detectorPosX");
    man->CreateNtupleDColumn("detectorPosY");
    man->CreateNtupleDColumn("detectorPosZ");
    man->CreateNtupleDColumn("globalTime");
    man->CreateNtupleDColumn("creationTime");
    man->FinishNtuple(4);


    man->CreateNtuple("Scattering", "Scattering");
    man->CreateNtupleDColumn("posX");
    man->CreateNtupleDColumn("posY");
    man->CreateNtupleDColumn("posZ");
    man->CreateNtupleDColumn("scatterTime");
    man->CreateNtupleDColumn("edep");
    man->CreateNtupleSColumn("procName");
    man->CreateNtupleIColumn("singleScatter");
    man->FinishNtuple(5);

    man->CreateNtuple("Electrons", "Electrons");
    man->CreateNtupleDColumn("posX");
    man->CreateNtupleDColumn("posY");
    man->CreateNtupleDColumn("posZ");
    man->CreateNtupleDColumn("scatterTime");
    man->CreateNtupleDColumn("edep");
    man->FinishNtuple(6);
}

MyRunAction::~MyRunAction() {}


void MyRunAction::BeginOfRunAction(const G4Run* run)


{
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    G4int runID = run->GetRunID();

    std::stringstream strRunID;
    strRunID << runID;

    man->OpenFile("output"+strRunID.str()+".root");

}

void MyRunAction::EndOfRunAction(const G4Run*) // <-- Added "MyRunAction::"
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->Write();
    man->CloseFile();
}
                                                          
