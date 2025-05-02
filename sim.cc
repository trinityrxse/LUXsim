#include <iostream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "construction.hh"
#include "DMXPhysicsList.hh"
#include "action.hh"

#include <fstream>  

int main(int argc, char** argv)
{
    // Redirect all output to null stream
    std::ofstream nullstream("/dev/null");  // On Linux/MacOS, "NUL" on Windows
    std::streambuf* coutbuf = std::cout.rdbuf();  // Save the original stream buffer
    std::streambuf* cerrbuf = std::cerr.rdbuf();  // Save the original error buffer

    std::cout.rdbuf(nullstream.rdbuf());  // Redirect G4cout to null stream
    std::cerr.rdbuf(nullstream.rdbuf());  // Redirect G4cerr to null stream

    G4RunManager *runManager = new G4RunManager();
    runManager->SetUserInitialization(new MyDetectorConstruction());
    runManager->SetUserInitialization(new DMXPhysicsList());
    runManager->SetUserInitialization(new MyActionInitialization());
    runManager->Initialize();

    G4UIExecutive *ui = 0;
    
    if(argc == 1)
    {
        ui = new G4UIExecutive(argc, argv);
    }

    G4VisManager *visManager = new G4VisExecutive();
    visManager->Initialize();

    G4UImanager *UImanager = G4UImanager::GetUIpointer();
    if(ui)
    {
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
    }
    else
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }

    //UImanager->ApplyCommand("/control/execute vis.mac");

    ui->SessionStart();
    // delete visManager;
    // delete runManager;

    std::cout.rdbuf(coutbuf);
    std::cerr.rdbuf(cerrbuf);

    return 0;
}