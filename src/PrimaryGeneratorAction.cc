#include "PrimaryGeneratorAction.hh"

#include "EventInformation.hh"

#include <G4Event.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <Randomize.hh>
#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction() {
  particle_gun_ = new G4ParticleGun(1);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete particle_gun_;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  const auto energy_gev = SampleEnergyGeV();
  const auto energy = energy_gev * GeV;

  const auto cos_theta = SampleCosTheta();
  const auto sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
  const auto phi = 2.0 * CLHEP::pi * G4UniformRand();

  const auto direction = G4ThreeVector(sin_theta * std::cos(phi),
                                       sin_theta * std::sin(phi),
                                       -cos_theta);

  const auto plane_half = 40.0 * cm;
  const auto x0 = (G4UniformRand() * 2.0 - 1.0) * plane_half;
  const auto y0 = (G4UniformRand() * 2.0 - 1.0) * plane_half;
  auto position = G4ThreeVector(x0, y0, 2.0 * m);

  const double mu_plus_ratio = 1.25;
  const bool is_mu_plus = (G4UniformRand() < (mu_plus_ratio / (1.0 + mu_plus_ratio)));
  auto particle_name = is_mu_plus ? "mu+" : "mu-";
  auto particle_def = G4ParticleTable::GetParticleTable()->FindParticle(particle_name);

  particle_gun_->SetParticleDefinition(particle_def);
  particle_gun_->SetParticleEnergy(energy);
  particle_gun_->SetParticleMomentumDirection(direction);
  particle_gun_->SetParticlePosition(position);

  auto info = new EventInformation();
  info->SetPrimaryEnergy(energy_gev);
  info->SetPrimaryDirection(direction);
  info->SetPrimaryPosition(position);
  info->SetCharge(is_mu_plus ? 1 : -1);
  event->SetUserInformation(info);

  particle_gun_->GeneratePrimaryVertex(event);
}

double PrimaryGeneratorAction::SampleEnergyGeV() const {
  const double e_min = 1.0;
  const double e_max = 100.0;
  const double gamma = 2.7;

  const double pow_min = std::pow(e_min, 1.0 - gamma);
  const double pow_max = std::pow(e_max, 1.0 - gamma);
  const double rnd = G4UniformRand();
  const double energy = std::pow(pow_min + rnd * (pow_max - pow_min), 1.0 / (1.0 - gamma));
  return energy;
}

double PrimaryGeneratorAction::SampleCosTheta() const {
  while (true) {
    const double cos_theta = G4UniformRand();
    const double accept = cos_theta * cos_theta;
    if (G4UniformRand() < accept) {
      return cos_theta;
    }
  }
}
