#include "RunAction.hh"

#include "AnalysisManager.hh"

#include <G4Run.hh>

void RunAction::BeginOfRunAction(const G4Run*) {
  auto analysis = AnalysisManager::Instance();
  analysis->Book();
  analysis->BeginOfRun();
}

void RunAction::EndOfRunAction(const G4Run*) {
  auto analysis = AnalysisManager::Instance();
  analysis->EndOfRun();
  AnalysisManager::Destroy();
}
