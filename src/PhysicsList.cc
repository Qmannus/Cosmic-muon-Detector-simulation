#include "PhysicsList.hh"

#include <G4DecayPhysics.hh>
#include <G4EmStandardPhysics.hh>
#include <G4HadronPhysicsFTFP_BERT.hh>
#include <G4IonPhysics.hh>
#include <G4SystemOfUnits.hh>

PhysicsList::PhysicsList() {
  SetVerboseLevel(1);
  RegisterPhysics(new G4EmStandardPhysics());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
  RegisterPhysics(new G4IonPhysics());
}

void PhysicsList::SetCuts() {
  SetCutsWithDefault();
}
