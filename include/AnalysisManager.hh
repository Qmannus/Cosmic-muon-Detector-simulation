#ifndef ANALYSIS_MANAGER_HH
#define ANALYSIS_MANAGER_HH

#include <G4ThreeVector.hh>

class G4AnalysisManager;

class AnalysisManager {
public:
  static AnalysisManager* Instance();
  static void Destroy();

  void Book();
  void Save();
  void BeginOfRun();
  void EndOfRun();

  void FillEvent(double energy_gev, double cos_theta, double phi,
                 int charge, bool detected, double time_diff_ns,
                 double reco_x_cm, double reco_y_cm, double chi2);

  void FillHit(double x_cm, double y_cm, double z_cm, double time_ns, double edep_mev);
  void FillPositionResidual(double residual_cm);
  void FillSiPM(int sipm_id, double photons, double time_ns);
  void FillTrack(double theta, double phi, double chi2, double residual_cm);
  void FillEnergyDeposit(double edep_mev);
  void FillPhotonYield(double photons);
  void FillSiPMMultiplicity(int multiplicity);

  void FillEfficiencyEnergy(double energy_gev, bool detected);
  void FillEfficiencyAngle(double cos_theta, bool detected);

private:
  AnalysisManager();
  ~AnalysisManager() = default;

  G4AnalysisManager* manager_ = nullptr;
};

#endif
