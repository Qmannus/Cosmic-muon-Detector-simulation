#ifndef SIPM_SD_HH
#define SIPM_SD_HH

#include <G4VSensitiveDetector.hh>

class SiPMSD : public G4VSensitiveDetector {
public:
  explicit SiPMSD(const G4String& name);
  ~SiPMSD() override = default;

  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
};

#endif
