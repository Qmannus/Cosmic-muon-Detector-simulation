#include "ScintillatorSD.hh"

#include "EventAction.hh"

#include <G4EventManager.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4TouchableHistory.hh>

ScintillatorSD::ScintillatorSD(const G4String& name) : G4VSensitiveDetector(name) {}

G4bool ScintillatorSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  const auto edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.0) {
    return false;
  }

  auto event_action = static_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction());
  if (!event_action) {
    return false;
  }

  const auto touchable = step->GetPreStepPoint()->GetTouchable();
  const auto layer = touchable->GetCopyNumber();
  const auto pos = step->GetPreStepPoint()->GetPosition();
  const auto time_ns = step->GetPreStepPoint()->GetGlobalTime() / ns;
  const auto edep_mev = edep / MeV;

  event_action->AddScintillatorHit(layer, pos, time_ns, edep_mev);
  return true;
}
