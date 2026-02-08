#ifndef SCINTILLATOR_SD_HH
#define SCINTILLATOR_SD_HH

#include <G4VSensitiveDetector.hh>

class ScintillatorSD : public G4VSensitiveDetector {
public:
  explicit ScintillatorSD(const G4String& name);
  ~ScintillatorSD() override = default;

  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
};

#endif
