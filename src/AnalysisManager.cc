#include "AnalysisManager.hh"

#include <G4AnalysisManager.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4Threading.hh>

namespace {
G4ThreadLocal AnalysisManager* g_instance = nullptr;
}

AnalysisManager* AnalysisManager::Instance() {
  if (!g_instance) {
    g_instance = new AnalysisManager();
  }
  return g_instance;
}

void AnalysisManager::Destroy() {
  delete g_instance;
  g_instance = nullptr;
}

AnalysisManager::AnalysisManager() {
  manager_ = G4AnalysisManager::Instance();
}

void AnalysisManager::Book() {
  manager_->SetVerboseLevel(1);
  manager_->SetDefaultFileType("root");
  manager_->SetFileName("cosmic_muon");

  manager_->CreateH1("energy_spectrum", "Primary energy (GeV)", 200, 0.0, 100.0);
  manager_->CreateH1("energy_deposit", "Energy deposition per event (MeV)", 200, 0.0, 200.0);
  manager_->CreateH1("photon_yield", "Photon yield proxy", 200, 0.0, 5000.0);
  manager_->CreateH1("cos_theta", "cos(theta)", 120, 0.0, 1.0);
  manager_->CreateH1("azimuth", "phi (rad)", 120, -CLHEP::pi, CLHEP::pi);
  manager_->CreateH2("hit_map", "Hit map (cm)", 100, -50, 50, 100, -50, 50);
  manager_->CreateH1("position_residual", "Hit position residual (cm)", 200, -10, 10);
  manager_->CreateH1("chi2", "Track chi2", 200, 0.0, 50.0);
  manager_->CreateH1("track_residual", "Track residual (cm)", 200, 0.0, 10.0);
  manager_->CreateH1("sipm_photons", "SiPM photon proxy", 200, 0.0, 2000.0);
  manager_->CreateH1("sipm_multiplicity", "SiPM multiplicity", 10, -0.5, 9.5);
  manager_->CreateH1("time_diff", "Timing difference (ns)", 200, -20.0, 20.0);
  manager_->CreateH1("reco_x", "Reconstructed X (cm)", 200, -50, 50);
  manager_->CreateH1("reco_y", "Reconstructed Y (cm)", 200, -50, 50);
  manager_->CreateH1("track_theta", "Track theta (rad)", 180, 0.0, CLHEP::pi / 2);
  manager_->CreateH1("track_phi", "Track phi (rad)", 180, -CLHEP::pi, CLHEP::pi);
  manager_->CreateH1("eff_energy_total", "Efficiency total vs energy (GeV)", 50, 0.0, 100.0);
  manager_->CreateH1("eff_energy_detected", "Efficiency detected vs energy (GeV)", 50, 0.0, 100.0);
  manager_->CreateH1("eff_angle_total", "Efficiency total vs cos(theta)", 50, 0.0, 1.0);
  manager_->CreateH1("eff_angle_detected", "Efficiency detected vs cos(theta)", 50, 0.0, 1.0);

  manager_->CreateNtuple("Events", "Event-level data");
  manager_->CreateNtupleDColumn("energy_gev");
  manager_->CreateNtupleDColumn("cos_theta");
  manager_->CreateNtupleDColumn("phi");
  manager_->CreateNtupleIColumn("charge");
  manager_->CreateNtupleIColumn("detected");
  manager_->CreateNtupleDColumn("time_diff_ns");
  manager_->CreateNtupleDColumn("reco_x_cm");
  manager_->CreateNtupleDColumn("reco_y_cm");
  manager_->CreateNtupleDColumn("chi2");
  manager_->FinishNtuple();

  manager_->CreateNtuple("Hits", "Hit positions");
  manager_->CreateNtupleDColumn("x_cm");
  manager_->CreateNtupleDColumn("y_cm");
  manager_->CreateNtupleDColumn("z_cm");
  manager_->CreateNtupleDColumn("time_ns");
  manager_->CreateNtupleDColumn("edep_mev");
  manager_->FinishNtuple();

  manager_->CreateNtuple("SiPMs", "SiPM responses");
  manager_->CreateNtupleIColumn("sipm_id");
  manager_->CreateNtupleDColumn("photons");
  manager_->CreateNtupleDColumn("time_ns");
  manager_->FinishNtuple();

  manager_->CreateNtuple("Tracks", "Track parameters");
  manager_->CreateNtupleDColumn("theta");
  manager_->CreateNtupleDColumn("phi");
  manager_->CreateNtupleDColumn("chi2");
  manager_->CreateNtupleDColumn("residual_cm");
  manager_->FinishNtuple();
}

void AnalysisManager::BeginOfRun() {
  manager_->OpenFile();
}

void AnalysisManager::EndOfRun() {
  Save();
}

void AnalysisManager::Save() {
  manager_->Write();
  manager_->CloseFile();
}

void AnalysisManager::FillEvent(double energy_gev, double cos_theta, double phi,
                                int charge, bool detected, double time_diff_ns,
                                double reco_x_cm, double reco_y_cm, double chi2) {
  manager_->FillH1(0, energy_gev);
  manager_->FillH1(3, cos_theta);
  manager_->FillH1(4, phi);
  manager_->FillH1(10, time_diff_ns);
  manager_->FillH1(11, reco_x_cm);
  manager_->FillH1(12, reco_y_cm);
  manager_->FillH1(6, chi2);

  manager_->FillNtupleDColumn(0, 0, energy_gev);
  manager_->FillNtupleDColumn(0, 1, cos_theta);
  manager_->FillNtupleDColumn(0, 2, phi);
  manager_->FillNtupleIColumn(0, 3, charge);
  manager_->FillNtupleIColumn(0, 4, detected ? 1 : 0);
  manager_->FillNtupleDColumn(0, 5, time_diff_ns);
  manager_->FillNtupleDColumn(0, 6, reco_x_cm);
  manager_->FillNtupleDColumn(0, 7, reco_y_cm);
  manager_->FillNtupleDColumn(0, 8, chi2);
  manager_->AddNtupleRow(0);
}

void AnalysisManager::FillHit(double x_cm, double y_cm, double z_cm, double time_ns, double edep_mev) {
  manager_->FillH2(0, x_cm, y_cm);

  manager_->FillNtupleDColumn(1, 0, x_cm);
  manager_->FillNtupleDColumn(1, 1, y_cm);
  manager_->FillNtupleDColumn(1, 2, z_cm);
  manager_->FillNtupleDColumn(1, 3, time_ns);
  manager_->FillNtupleDColumn(1, 4, edep_mev);
  manager_->AddNtupleRow(1);
}

void AnalysisManager::FillSiPM(int sipm_id, double photons, double time_ns) {
  manager_->FillH1(8, photons);

  manager_->FillNtupleIColumn(2, 0, sipm_id);
  manager_->FillNtupleDColumn(2, 1, photons);
  manager_->FillNtupleDColumn(2, 2, time_ns);
  manager_->AddNtupleRow(2);
}

void AnalysisManager::FillPositionResidual(double residual_cm) {
  manager_->FillH1(5, residual_cm);
}

void AnalysisManager::FillTrack(double theta, double phi, double chi2, double residual_cm) {
  manager_->FillH1(13, theta);
  manager_->FillH1(14, phi);
  manager_->FillH1(6, chi2);
  manager_->FillH1(7, residual_cm);

  manager_->FillNtupleDColumn(3, 0, theta);
  manager_->FillNtupleDColumn(3, 1, phi);
  manager_->FillNtupleDColumn(3, 2, chi2);
  manager_->FillNtupleDColumn(3, 3, residual_cm);
  manager_->AddNtupleRow(3);
}

void AnalysisManager::FillEnergyDeposit(double edep_mev) {
  manager_->FillH1(1, edep_mev);
}

void AnalysisManager::FillPhotonYield(double photons) {
  manager_->FillH1(2, photons);
}

void AnalysisManager::FillSiPMMultiplicity(int multiplicity) {
  manager_->FillH1(9, multiplicity);
}

void AnalysisManager::FillEfficiencyEnergy(double energy_gev, bool detected) {
  manager_->FillH1(15, energy_gev);
  if (detected) {
    manager_->FillH1(16, energy_gev);
  }
}

void AnalysisManager::FillEfficiencyAngle(double cos_theta, bool detected) {
  manager_->FillH1(17, cos_theta);
  if (detected) {
    manager_->FillH1(18, cos_theta);
  }
}
