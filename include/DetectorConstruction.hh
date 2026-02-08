#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include <G4VUserDetectorConstruction.hh>
#include <G4ThreeVector.hh>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction() override = default;

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  const std::vector<G4ThreeVector>& GetScintillatorCenters() const;
  const std::vector<G4ThreeVector>& GetSiPMCenters() const;

private:
  void DefineMaterials();

  G4Material* air_ = nullptr;
  G4Material* scintillator_ = nullptr;
  G4Material* silicon_ = nullptr;

  G4LogicalVolume* scintillator_logical_ = nullptr;
  G4LogicalVolume* sipm_logical_ = nullptr;

  std::vector<G4ThreeVector> scintillator_centers_;
  std::vector<G4ThreeVector> sipm_centers_;
};

#endif
