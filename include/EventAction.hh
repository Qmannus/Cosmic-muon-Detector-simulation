#ifndef EVENT_ACTION_HH
#define EVENT_ACTION_HH

#include <G4UserEventAction.hh>
#include <G4ThreeVector.hh>
#include <vector>

struct ScintHitData {
  int layer = 0;
  G4ThreeVector position;
  double time_ns = 0.0;
  double edep_mev = 0.0;
};

struct SiPMHitData {
  int sipm_id = 0;
  double photons = 0.0;
  double time_ns = 0.0;
};

class EventAction : public G4UserEventAction {
public:
  EventAction() = default;
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;

  void AddScintillatorHit(int layer, const G4ThreeVector& pos, double time_ns, double edep_mev);
  void AddSiPMHit(int sipm_id, double photons, double time_ns);

private:
  std::vector<ScintHitData> scint_hits_;
  std::vector<SiPMHitData> sipm_hits_;
  double total_edep_mev_ = 0.0;
  double total_photons_ = 0.0;
};

#endif
