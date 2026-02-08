#include "DetectorConstruction.hh"

#include "ScintillatorSD.hh"
#include "SiPMSD.hh"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VisAttributes.hh>

DetectorConstruction::DetectorConstruction() = default;

void DetectorConstruction::DefineMaterials() {
  auto nist = G4NistManager::Instance();
  air_ = nist->FindOrBuildMaterial("G4_AIR");
  scintillator_ = nist->FindOrBuildMaterial("G4_POLYVINYL_TOLUENE");
  if (!scintillator_) {
    scintillator_ = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  }
  silicon_ = nist->FindOrBuildMaterial("G4_Si");
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  DefineMaterials();

  auto world_half = 2.5 * m;
  auto world_solid = new G4Box("World", world_half, world_half, world_half);
  auto world_logical = new G4LogicalVolume(world_solid, air_, "WorldLogical");
  auto world_physical = new G4PVPlacement(nullptr, {}, world_logical, "World", nullptr, false, 0, true);

  auto scintillator_xy = 31.5 * cm;
  auto scintillator_z = 0.5 * cm;
  auto scintillator_solid = new G4Box("Scintillator", scintillator_xy, scintillator_xy, scintillator_z);
  scintillator_logical_ = new G4LogicalVolume(scintillator_solid, scintillator_, "ScintillatorLogical");

  const G4double layer_gap = 63.0 * cm;
  scintillator_centers_.clear();
  scintillator_centers_.push_back({0, 0, 0});
  scintillator_centers_.push_back({0, 0, layer_gap});

  for (size_t i = 0; i < scintillator_centers_.size(); ++i) {
    new G4PVPlacement(nullptr, scintillator_centers_[i], scintillator_logical_,
                      "Scintillator", world_logical, false, static_cast<int>(i), true);
  }

  auto sipm_half = 0.3 * cm;
  auto sipm_solid = new G4Box("SiPM", sipm_half, sipm_half, sipm_half);
  sipm_logical_ = new G4LogicalVolume(sipm_solid, silicon_, "SiPMLogical");

  sipm_centers_.clear();
  const G4double sipm_offset = scintillator_xy - sipm_half;
  std::vector<G4ThreeVector> sipm_offsets = {
      {sipm_offset, sipm_offset, 0},
      {-sipm_offset, sipm_offset, 0},
      {sipm_offset, -sipm_offset, 0},
      {-sipm_offset, -sipm_offset, 0},
      {sipm_offset, sipm_offset, layer_gap},
      {-sipm_offset, sipm_offset, layer_gap},
      {sipm_offset, -sipm_offset, layer_gap},
      {-sipm_offset, -sipm_offset, layer_gap}};

  for (size_t i = 0; i < sipm_offsets.size(); ++i) {
    sipm_centers_.push_back(sipm_offsets[i]);
    new G4PVPlacement(nullptr, sipm_offsets[i], sipm_logical_, "SiPM",
                      world_logical, false, static_cast<int>(i), true);
  }

  auto scint_attr = new G4VisAttributes(G4Colour(0.1, 0.6, 1.0, 0.35));
  scint_attr->SetForceSolid(true);
  scintillator_logical_->SetVisAttributes(scint_attr);

  auto sipm_attr = new G4VisAttributes(G4Colour(1.0, 0.1, 0.1));
  sipm_attr->SetForceSolid(true);
  sipm_logical_->SetVisAttributes(sipm_attr);

  auto world_attr = new G4VisAttributes();
  world_attr->SetVisibility(false);
  world_logical->SetVisAttributes(world_attr);

  return world_physical;
}

void DetectorConstruction::ConstructSDandField() {
  auto sd_manager = G4SDManager::GetSDMpointer();

  auto scint_sd = new ScintillatorSD("ScintillatorSD");
  sd_manager->AddNewDetector(scint_sd);
  SetSensitiveDetector(scintillator_logical_, scint_sd);

  auto sipm_sd = new SiPMSD("SiPMSD");
  sd_manager->AddNewDetector(sipm_sd);
  SetSensitiveDetector(sipm_logical_, sipm_sd);
}

const std::vector<G4ThreeVector>& DetectorConstruction::GetScintillatorCenters() const {
  return scintillator_centers_;
}

const std::vector<G4ThreeVector>& DetectorConstruction::GetSiPMCenters() const {
  return sipm_centers_;
}
