#include "EventAction.hh"

#include "AnalysisManager.hh"
#include "EventInformation.hh"

#include <G4Event.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <algorithm>
#include <cmath>
#include <limits>

namespace {
constexpr double kScintillatorHalfSizeCm = 31.5;
constexpr double kPlaneSeparationCm = 63.0;
constexpr double kRefractiveIndex = 1.58;
constexpr double kLightYieldPerMeV = 10000.0;
constexpr double kPDE = 0.41;
constexpr double kTriggerThresholdPE = 10.0;
constexpr double kSigmaElectronicsNs = 0.3354;

struct PlaneReconstruction {
  bool valid = false;
  G4ThreeVector poi = {};
  double t0_ns = 0.0;
  double chi2 = 0.0;
  std::vector<SiPMHitData> sipm_hits;
};

std::vector<G4ThreeVector> MakePlaneSiPMPositions(double z_cm) {
  const double offset = kScintillatorHalfSizeCm - 0.3;
  return {
      {offset * cm, offset * cm, z_cm * cm},
      {-offset * cm, offset * cm, z_cm * cm},
      {offset * cm, -offset * cm, z_cm * cm},
      {-offset * cm, -offset * cm, z_cm * cm},
  };
}

PlaneReconstruction ReconstructPlane(const G4ThreeVector& hit_pos, double hit_time_ns,
                                     double plane_z_cm, double edep_mev, int sipm_id_offset) {
  PlaneReconstruction result;
  const auto sipm_positions = MakePlaneSiPMPositions(plane_z_cm);

  const double photons = edep_mev * kLightYieldPerMeV;
  const double expected_pe = photons * kPDE / 4.0;
  const double v_gamma = (CLHEP::c_light / kRefractiveIndex) / (cm / ns);

  std::vector<double> measured_times;
  measured_times.reserve(4);

  for (size_t i = 0; i < sipm_positions.size(); ++i) {
    const auto distance_cm = (hit_pos - sipm_positions[i]).mag() / cm;
    const double tof_ns = distance_cm / v_gamma;
    const double pe = G4Poisson(expected_pe);
    if (pe < kTriggerThresholdPE) {
      measured_times.push_back(std::numeric_limits<double>::quiet_NaN());
      continue;
    }
    const double smeared_time = hit_time_ns + tof_ns + G4RandGauss::shoot(0.0, kSigmaElectronicsNs);
    measured_times.push_back(smeared_time);
    result.sipm_hits.push_back({static_cast<int>(sipm_id_offset + i), pe, smeared_time});
  }

  int valid_hits = 0;
  for (double t : measured_times) {
    if (std::isfinite(t)) {
      valid_hits++;
    }
  }
  if (valid_hits < 3) {
    return result;
  }

  double best_chi2 = 1e30;
  G4ThreeVector best_poi;
  double best_t0 = 0.0;

  for (double x_cm = -kScintillatorHalfSizeCm; x_cm <= kScintillatorHalfSizeCm; x_cm += 1.0) {
    for (double y_cm = -kScintillatorHalfSizeCm; y_cm <= kScintillatorHalfSizeCm; y_cm += 1.0) {
      const auto poi = G4ThreeVector(x_cm * cm, y_cm * cm, plane_z_cm * cm);
      double t0_sum = 0.0;
      int t0_count = 0;
      for (size_t i = 0; i < sipm_positions.size(); ++i) {
        if (!std::isfinite(measured_times[i])) {
          continue;
        }
        const auto distance_cm = (poi - sipm_positions[i]).mag() / cm;
        const double tof_ns = distance_cm / v_gamma;
        t0_sum += measured_times[i] - tof_ns;
        t0_count++;
      }
      if (t0_count == 0) {
        continue;
      }
      const double t0 = t0_sum / t0_count;
      double chi2 = 0.0;
      for (size_t i = 0; i < sipm_positions.size(); ++i) {
        if (!std::isfinite(measured_times[i])) {
          continue;
        }
        const auto distance_cm = (poi - sipm_positions[i]).mag() / cm;
        const double tof_ns = distance_cm / v_gamma;
        const double residual = measured_times[i] - (t0 + tof_ns);
        chi2 += (residual * residual) / (kSigmaElectronicsNs * kSigmaElectronicsNs);
      }
      if (chi2 < best_chi2) {
        best_chi2 = chi2;
        best_poi = poi;
        best_t0 = t0;
      }
    }
  }

  result.valid = true;
  result.poi = best_poi;
  result.t0_ns = best_t0;
  result.chi2 = best_chi2;
  return result;
}
}  // namespace

void EventAction::BeginOfEventAction(const G4Event*) {
  scint_hits_.clear();
  sipm_hits_.clear();
  total_edep_mev_ = 0.0;
  total_photons_ = 0.0;
}

void EventAction::EndOfEventAction(const G4Event* event) {
  auto analysis = AnalysisManager::Instance();

  bool has_layer0 = false;
  bool has_layer1 = false;
  double earliest_layer0 = 1e9;
  double earliest_layer1 = 1e9;
  G4ThreeVector sum_layer0;
  G4ThreeVector sum_layer1;
  int count_layer0 = 0;
  int count_layer1 = 0;
  double edep_layer0 = 0.0;
  double edep_layer1 = 0.0;

  for (const auto& hit : scint_hits_) {
    analysis->FillHit(hit.position.x() / cm, hit.position.y() / cm, hit.position.z() / cm,
                      hit.time_ns, hit.edep_mev);

    if (hit.layer == 0) {
      has_layer0 = true;
      earliest_layer0 = std::min(earliest_layer0, hit.time_ns);
      sum_layer0 += hit.position;
      count_layer0++;
      edep_layer0 += hit.edep_mev;
    } else {
      has_layer1 = true;
      earliest_layer1 = std::min(earliest_layer1, hit.time_ns);
      sum_layer1 += hit.position;
      count_layer1++;
      edep_layer1 += hit.edep_mev;
    }
  }

  PlaneReconstruction plane0;
  PlaneReconstruction plane1;
  if (has_layer0 && count_layer0 > 0) {
    auto avg0 = sum_layer0 / count_layer0;
    plane0 = ReconstructPlane(avg0, earliest_layer0, 0.0, edep_layer0, 0);
  }
  if (has_layer1 && count_layer1 > 0) {
    auto avg1 = sum_layer1 / count_layer1;
    plane1 = ReconstructPlane(avg1, earliest_layer1, kPlaneSeparationCm, edep_layer1, 4);
  }

  total_photons_ = (edep_layer0 + edep_layer1) * kLightYieldPerMeV;
  analysis->FillPhotonYield(total_photons_);
  analysis->FillSiPMMultiplicity(static_cast<int>(plane0.sipm_hits.size() + plane1.sipm_hits.size()));
  for (const auto& hit : plane0.sipm_hits) {
    analysis->FillSiPM(hit.sipm_id, hit.photons, hit.time_ns);
  }
  for (const auto& hit : plane1.sipm_hits) {
    analysis->FillSiPM(hit.sipm_id, hit.photons, hit.time_ns);
  }

  bool detected = plane0.valid && plane1.valid;
  double time_diff = 0.0;
  if (detected) {
    time_diff = plane1.t0_ns - plane0.t0_ns;
  }

  G4ThreeVector reco_point;
  if (detected) {
    reco_point = (plane0.poi + plane1.poi) * 0.5;
  }

  double chi2 = 0.0;
  double residual_cm = 0.0;
  double theta = 0.0;
  double phi = 0.0;

  const auto* info = dynamic_cast<EventInformation*>(event->GetUserInformation());
  double energy_gev = info ? info->GetPrimaryEnergy() : 0.0;
  int charge = info ? info->GetCharge() : 0;
  auto true_dir = info ? info->GetPrimaryDirection() : G4ThreeVector(0, 0, -1);
  auto true_pos = info ? info->GetPrimaryPosition() : G4ThreeVector();
  double cos_theta = std::abs(true_dir.z());
  double true_phi = std::atan2(true_dir.y(), true_dir.x());

  if (detected) {
    auto reco_dir = (plane1.poi - plane0.poi).unit();
    theta = std::acos(std::abs(reco_dir.z()));
    phi = std::atan2(reco_dir.y(), reco_dir.x());

    if (std::abs(true_dir.z()) > 1e-6) {
      double t = (0.0 - true_pos.z()) / true_dir.z();
      auto true_mid = true_pos + t * true_dir;
      residual_cm = (reco_point - true_mid).mag() / cm;
    }

    chi2 = plane0.chi2 + plane1.chi2;

    analysis->FillTrack(theta, phi, chi2, residual_cm);
    analysis->FillEfficiencyEnergy(energy_gev, true);
    analysis->FillEfficiencyAngle(cos_theta, true);

    if (std::abs(true_dir.z()) > 1e-6) {
      double t0 = (0.0 - true_pos.z()) / true_dir.z();
      auto true_plane0 = true_pos + t0 * true_dir;
      double t1 = (kPlaneSeparationCm * cm - true_pos.z()) / true_dir.z();
      auto true_plane1 = true_pos + t1 * true_dir;
      analysis->FillPositionResidual((plane0.poi - true_plane0).mag() / cm);
      analysis->FillPositionResidual((plane1.poi - true_plane1).mag() / cm);
    }
  } else {
    analysis->FillEfficiencyEnergy(energy_gev, false);
    analysis->FillEfficiencyAngle(cos_theta, false);
  }

  analysis->FillEvent(energy_gev, cos_theta, true_phi, charge, detected, time_diff,
                      reco_point.x() / cm, reco_point.y() / cm, chi2);

  analysis->FillEnergyDeposit(total_edep_mev_);
}

void EventAction::AddScintillatorHit(int layer, const G4ThreeVector& pos, double time_ns, double edep_mev) {
  scint_hits_.push_back({layer, pos, time_ns, edep_mev});
  total_edep_mev_ += edep_mev;
}

void EventAction::AddSiPMHit(int sipm_id, double photons, double time_ns) {
  sipm_hits_.push_back({sipm_id, photons, time_ns});
  total_photons_ += photons;
}
