# Cosmic Muon Detector Simulation (Geant4 11.4.0)

This repository provides a Geant4-based simulation of a two-layer scintillator cosmic muon detector with SiPM readout, ROOT output, and a Python analysis script that generates publication-ready plots.

## Features
- Cosmic muon primary generator with $E^{-2.7}$ spectrum and $\cos^2\theta$ angular distribution.
- Semi-transparent 63x63x1 cm scintillator panels (BC-408-like) separated by 63 cm, with highlighted 6 mm SiPMs.
- 3D visualization macro that accumulates 50 events.
- Electronics smearing with 41% PDE, 10 PE trigger threshold, and 335.4 ps timing resolution.
- ROOT output with 20 histograms and 4 ntuples (Events, Hits, SiPMs, Tracks).
- Python analysis script to generate 7 summary plots.

## Build (CMake + Visual Studio)
```bash
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022" -A x64
cmake --build . --config Release
cmake --build . --config Release --target INSTALL
```

## Run
Interactive visualization:
```bash
./cosmic_muon_sim
```

Batch mode (writes `cosmic_muon.root`):
```bash
./cosmic_muon_sim macros/run.mac
```

## Analysis
```bash
python python/analyze_data.py
```

Plots are saved to `analysis_plots/`.
