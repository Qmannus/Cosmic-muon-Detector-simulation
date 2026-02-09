#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include <G4ios.hh>
#include <G4RunManagerFactory.hh>
#include <G4UImanager.hh>

int main(int argc, char** argv) {
  auto run_manager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  run_manager->SetUserInitialization(new DetectorConstruction());
  run_manager->SetUserInitialization(new PhysicsList());
  run_manager->SetUserInitialization(new ActionInitialization());

  run_manager->Initialize();

  auto ui_manager = G4UImanager::GetUIpointer();
  ui_manager->ApplyCommand("/control/macroPath macros");
  if (argc > 1) {
    auto command = G4String("/control/execute ");
    G4String file_name = argv[1];
    ui_manager->ApplyCommand(command + file_name);
  } else {
    ui_manager->ApplyCommand("/control/execute run.mac");
  }

  delete run_manager;
  return 0;
}
