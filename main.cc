#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Shielding.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "EventAction.hh"

int main(int argc, char** argv) {
    // Detect interactive mode (if no arguments) and define UI session
    G4UIExecutive* ui = nullptr;
    if (argc == 1 && std::getenv("BATCH_MODE") == nullptr) {
        ui = new G4UIExecutive(argc, argv); // GUI mode only if not in batch
    }

    // Construct the default run manager
    G4RunManager* runManager = new G4RunManager;

    // Set mandatory initialization classes
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new Shielding());

    // Initialize the visualization (only needed for GUI mode)
    //G4VisManager* visManager = nullptr;
    //if (ui) {
    //    visManager = new G4VisExecutive;
    //    visManager->Initialize();
    //}
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    // Set user action classes
    runManager->SetUserAction(new PrimaryGeneratorAction());


     // Create RunAction and EventAction
    RunAction* runAction = new RunAction();
    EventAction* eventAction = new EventAction();

    // Connect RunAction to EventAction
    eventAction->SetRunAction(runAction);

    // Set the actions
    runManager->SetUserAction(runAction);
    runManager->SetUserAction(eventAction);

    // Pass the output directory to SteppingAction

    // for running from the project root 'muon-induced-secondaries'
    G4String outputDir = "./secondaries/secondary_data";

    // for running from the release directory with the .exe
    //G4String outputDir = "../../secondaries/secondary_data";

    if (argc > 2) {
        outputDir = argv[2]; // Output directory passed as the second argument
    }
    runManager->SetUserAction(new SteppingAction(outputDir));

    // Initialize G4 kernel
    runManager->Initialize();

    // Get the pointer to the UI manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (ui) {
        // Interactive mode (GUI)
        // uncomment to startup with macro command
        //UImanager->ApplyCommand("/control/execute C:/Users/User/GEANT4/share/Geant4/examples/muon-induced-secondaries/macros/init_vis.mac");
        ui->SessionStart();
        delete ui;
    }
    else {
        // Batch mode (no GUI)
        G4String command = "/control/execute ";
        G4String fileName = argv[1]; // Macro file passed as the first argument
        UImanager->ApplyCommand(command + fileName);
    }

    // Job termination
    if (visManager) delete visManager;
    delete runManager;
    return 0;
}