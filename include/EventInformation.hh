#ifndef EVENT_INFORMATION_HH
#define EVENT_INFORMATION_HH

#include <G4VUserEventInformation.hh>
#include <G4ThreeVector.hh>

class EventInformation : public G4VUserEventInformation {
public:
  void Print() const override {}

  void SetPrimaryEnergy(double energy) { primary_energy_ = energy; }
  void SetPrimaryDirection(const G4ThreeVector& dir) { primary_direction_ = dir; }
  void SetPrimaryPosition(const G4ThreeVector& pos) { primary_position_ = pos; }
  void SetCharge(int charge) { charge_ = charge; }

  double GetPrimaryEnergy() const { return primary_energy_; }
  const G4ThreeVector& GetPrimaryDirection() const { return primary_direction_; }
  const G4ThreeVector& GetPrimaryPosition() const { return primary_position_; }
  int GetCharge() const { return charge_; }

private:
  double primary_energy_ = 0.0;
  G4ThreeVector primary_direction_ = {0, 0, -1};
  G4ThreeVector primary_position_ = {};
  int charge_ = 0;
};

#endif
