cmake_minimum_required(VERSION 3.16...3.27)
project(CosmicMuonDetectorSimulation LANGUAGES CXX)

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and visualization" ON)

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/include)

set(SOURCES
  src/main.cc
  src/DetectorConstruction.cc
  src/PhysicsList.cc
  src/PrimaryGeneratorAction.cc
  src/RunAction.cc
  src/EventAction.cc
  src/ScintillatorSD.cc
  src/SiPMSD.cc
  src/AnalysisManager.cc
)

add_executable(cosmic_muon_sim ${SOURCES})

target_link_libraries(cosmic_muon_sim ${Geant4_LIBRARIES})

install(TARGETS cosmic_muon_sim DESTINATION bin)

file(GLOB MACROS "${PROJECT_SOURCE_DIR}/macros/*.mac")
install(FILES ${MACROS} DESTINATION bin)
install(FILES ${PROJECT_SOURCE_DIR}/python/analyze_data.py DESTINATION bin)
