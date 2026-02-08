#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include <G4RunManagerFactory.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>

int main(int argc, char** argv) {
  auto run_manager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  run_manager->SetUserInitialization(new DetectorConstruction());
  run_manager->SetUserInitialization(new PhysicsList());
  run_manager->SetUserInitialization(new ActionInitialization());

  run_manager->Initialize();

  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  auto vis_manager = new G4VisExecutive();
  vis_manager->Initialize();

  auto ui_manager = G4UImanager::GetUIpointer();
  if (ui) {
    ui_manager->ApplyCommand("/control/execute macros/vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    auto command = G4String("/control/execute ");
    G4String file_name = argv[1];
    ui_manager->ApplyCommand(command + file_name);
  }

  delete vis_manager;
  delete run_manager;
  return 0;
}
