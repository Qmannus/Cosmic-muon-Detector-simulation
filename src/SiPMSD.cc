#include "SiPMSD.hh"

#include <G4Step.hh>

SiPMSD::SiPMSD(const G4String& name) : G4VSensitiveDetector(name) {}

G4bool SiPMSD::ProcessHits(G4Step*, G4TouchableHistory*) {
  return false;
}
