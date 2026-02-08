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

  auto element_h = nist->FindOrBuildElement("H");
  auto element_c = nist->FindOrBuildElement("C");
  auto element_n = nist->FindOrBuildElement("N");
  auto element_o = nist->FindOrBuildElement("O");
  auto element_ar = nist->FindOrBuildElement("Ar");
  auto element_si = nist->FindOrBuildElement("Si");

  air_ = new G4Material("Air", 1.225 * mg / cm3, 3);
  air_->AddElement(element_n, 0.78);
  air_->AddElement(element_o, 0.21);
  air_->AddElement(element_ar, 0.01);

  scintillator_ = new G4Material("BC408", 1.023 * g / cm3, 2);
  scintillator_->AddElement(element_c, 1000);
  scintillator_->AddElement(element_h, 1104);

  silicon_ = new G4Material("Silicon", 2.33 * g / cm3, 1);
  silicon_->AddElement(element_si, 1);
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
