#!/usr/bin/env python3
import pathlib

import matplotlib.pyplot as plt
import numpy as np

try:
    import uproot
except ImportError as exc:
    raise SystemExit("This script requires uproot. Install with: pip install uproot") from exc

ROOT_FILE = pathlib.Path("cosmic_muon.root")
OUTPUT_DIR = pathlib.Path("analysis_plots")
OUTPUT_DIR.mkdir(exist_ok=True)

PLANE_SEPARATION_CM = 63.0


def _gaussian_fit(data):
    if len(data) == 0:
        return 0.0, 1.0
    mean = np.mean(data)
    sigma = np.std(data)
    return mean, sigma


def _hist_stats(hist):
    counts, edges = hist
    centers = 0.5 * (edges[1:] + edges[:-1])
    total = np.sum(counts)
    if total <= 0:
        return 0.0, 1.0
    mean = np.sum(centers * counts) / total
    variance = np.sum(counts * (centers - mean) ** 2) / total
    return mean, np.sqrt(variance)


def main():
    if not ROOT_FILE.exists():
        raise SystemExit(f"ROOT file '{ROOT_FILE}' not found. Run the simulation first.")

    with uproot.open(ROOT_FILE) as file:
        events = file["Events"].arrays(library="np")
        hits = file["Hits"].arrays(library="np")
        sipms = file["SiPMs"].arrays(library="np")
        tracks = file["Tracks"].arrays(library="np")
        h_energy = file["energy_spectrum"].to_numpy()
        h_edep = file["energy_deposit"].to_numpy()
        h_photon = file["photon_yield"].to_numpy()
        h_cos = file["cos_theta"].to_numpy()
        h_phi = file["azimuth"].to_numpy()
        h_hit_map = file["hit_map"].to_numpy()
        h_chi2 = file["chi2"].to_numpy()
        h_resid = file["position_residual"].to_numpy()
        h_sipm = file["sipm_photons"].to_numpy()
        h_eff_energy_total = file["eff_energy_total"].to_numpy()
        h_eff_energy_det = file["eff_energy_detected"].to_numpy()
        h_eff_angle_total = file["eff_angle_total"].to_numpy()
        h_eff_angle_det = file["eff_angle_detected"].to_numpy()

    # Energy distributions
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].stairs(h_energy[0], h_energy[1], fill=True, alpha=0.7)
    ax[0].set_xlabel("Energy [GeV]")
    ax[0].set_ylabel("Counts")
    ax[0].set_title("Primary spectrum")

    ax[1].stairs(h_edep[0], h_edep[1], fill=True, alpha=0.7, color="orange")
    ax[1].set_xlabel("Deposited energy [MeV]")
    ax[1].set_title("Total deposited energy")

    ax[2].stairs(h_photon[0], h_photon[1], fill=True, alpha=0.7, color="green")
    ax[2].set_xlabel("Photon yield proxy")
    ax[2].set_title("Photon yield")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "energy_distributions.png", dpi=200)
    plt.close(fig)

    # Angular distributions
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].stairs(h_cos[0], h_cos[1], fill=True, alpha=0.7)
    cos_centers = 0.5 * (h_cos[1][1:] + h_cos[1][:-1])
    ax[0].plot(cos_centers, np.max(h_cos[0]) * cos_centers**2, "k--", label=r"$\cos^2\theta$")
    ax[0].set_xlabel("cos(theta)")
    ax[0].set_ylabel("Counts")
    ax[0].legend()

    ax[1].stairs(h_phi[0], h_phi[1], fill=True, alpha=0.7, color="purple")
    ax[1].set_xlabel("phi [rad]")
    ax[1].set_ylabel("Counts")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "angular_distributions.png", dpi=200)
    plt.close(fig)

    # 2D hit map
    fig, ax = plt.subplots(figsize=(6, 6))
    x_edges, y_edges = h_hit_map[1], h_hit_map[2]
    ax.pcolormesh(x_edges, y_edges, h_hit_map[0].T, cmap="viridis")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.set_title("Hit map")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "2d_hitmap.png", dpi=200)
    plt.close(fig)

    # Reconstruction quality
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].stairs(h_chi2[0], h_chi2[1], fill=True, alpha=0.7)
    ax[0].set_xlabel("chi2")
    ax[0].set_ylabel("Counts")

    track_residuals = tracks["residual_cm"] if len(tracks) else np.array([])
    mean, sigma = _gaussian_fit(track_residuals)
    ax[1].hist(track_residuals, bins=60, color="slateblue", alpha=0.7, density=True)
    x_vals = np.linspace(mean - 4 * sigma, mean + 4 * sigma, 200)
    if sigma > 0:
        gauss = 1.0 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x_vals - mean) / sigma) ** 2)
        ax[1].plot(x_vals, gauss, "k--", label=f"track σ={sigma:.2f} cm")
    ax[1].set_xlabel("Track residual [cm]")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "reconstruction_quality.png", dpi=200)
    plt.close(fig)

    # SiPM performance
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].stairs(h_sipm[0], h_sipm[1], fill=True, alpha=0.7, color="red")
    ax[0].set_xlabel("Photon proxy")
    ax[0].set_ylabel("Counts")

    ax[1].hist(sipms["sipm_id"], bins=np.arange(-0.5, 8.5, 1), color="gray")
    ax[1].set_xlabel("SiPM ID")
    ax[1].set_ylabel("Multiplicity")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "sipm_performance.png", dpi=200)
    plt.close(fig)

    # Efficiency studies
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    eff_energy = np.divide(h_eff_energy_det[0], h_eff_energy_total[0],
                           out=np.zeros_like(h_eff_energy_det[0]), where=h_eff_energy_total[0] > 0)
    eff_angle = np.divide(h_eff_angle_det[0], h_eff_angle_total[0],
                          out=np.zeros_like(h_eff_angle_det[0]), where=h_eff_angle_total[0] > 0)
    energy_centers = 0.5 * (h_eff_energy_total[1][1:] + h_eff_energy_total[1][:-1])
    angle_centers = 0.5 * (h_eff_angle_total[1][1:] + h_eff_angle_total[1][:-1])
    ax[0].plot(energy_centers, eff_energy, "o-")
    ax[0].set_xlabel("Energy [GeV]")
    ax[0].set_ylabel("Efficiency")
    ax[1].plot(angle_centers, eff_angle, "o-")
    ax[1].set_xlabel("cos(theta)")
    ax[1].set_ylabel("Efficiency")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "efficiency_studies.png", dpi=200)
    plt.close(fig)

    # Timing + reconstruction distributions
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))
    ax = ax.flatten()
    h_time_diff = file["time_diff"].to_numpy()
    h_reco_x = file["reco_x"].to_numpy()
    h_reco_y = file["reco_y"].to_numpy()
    h_track_theta = file["track_theta"].to_numpy()

    ax[0].stairs(h_time_diff[0], h_time_diff[1], fill=True, alpha=0.7, color="teal")
    ax[0].set_xlabel("Δt [ns]")
    ax[0].set_ylabel("Counts")
    ax[0].set_title("Plane timing difference")

    ax[1].stairs(h_reco_x[0], h_reco_x[1], fill=True, alpha=0.7, color="steelblue")
    ax[1].set_xlabel("Reco X [cm]")
    ax[1].set_ylabel("Counts")
    ax[1].set_title("Reconstructed X")

    ax[2].stairs(h_reco_y[0], h_reco_y[1], fill=True, alpha=0.7, color="slateblue")
    ax[2].set_xlabel("Reco Y [cm]")
    ax[2].set_ylabel("Counts")
    ax[2].set_title("Reconstructed Y")

    ax[3].stairs(h_track_theta[0], h_track_theta[1], fill=True, alpha=0.7, color="coral")
    ax[3].set_xlabel("Track θ [rad]")
    ax[3].set_ylabel("Counts")
    ax[3].set_title("Track zenith angle")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "timing_reco_distributions.png", dpi=200)
    plt.close(fig)

    # Chi2 vs residuals
    fig, ax = plt.subplots(figsize=(6, 5))
    if len(tracks):
        ax.scatter(tracks["chi2"], tracks["residual_cm"], s=6, alpha=0.4)
    ax.set_xlabel("Track chi2")
    ax.set_ylabel("Residual [cm]")
    ax.set_title("Reconstruction quality correlation")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "chi2_vs_residual.png", dpi=200)
    plt.close(fig)

    # 3D reconstructed trajectories (using track angles + reconstructed midpoint)
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title("Reconstructed track samples")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.set_zlabel("z [cm]")

    if len(events) and len(tracks):
        z0 = 0.0
        z1 = PLANE_SEPARATION_CM
        z_mid = 0.5 * (z0 + z1)
        max_tracks = min(len(tracks), len(events))
        sample_count = min(max_tracks, 200)
        if sample_count > 0:
            indices = np.random.choice(max_tracks, size=sample_count, replace=False)
            for idx in indices:
                theta = tracks["theta"][idx]
                phi = tracks["phi"][idx]
                x_mid = events["reco_x_cm"][idx]
                y_mid = events["reco_y_cm"][idx]
                values = np.array([theta, phi, x_mid, y_mid], dtype=float)
                if not np.isfinite(values).all():
                    continue
                dx = np.sin(theta) * np.cos(phi)
                dy = np.sin(theta) * np.sin(phi)
                dz = np.cos(theta)
                if abs(dz) < 1e-6:
                    continue
                t0 = (z0 - z_mid) / dz
                t1 = (z1 - z_mid) / dz
                x0 = x_mid + dx * t0
                y0 = y_mid + dy * t0
                x1 = x_mid + dx * t1
                y1 = y_mid + dy * t1
                ax.plot([x0, x1], [y0, y1], [z0, z1], alpha=0.3, linewidth=0.7)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "trajectory_reconstruction.png", dpi=200)
    plt.close(fig)

    # Summary plot
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    axes = axes.flatten()
    axes[0].stairs(h_energy[0], h_energy[1])
    axes[0].set_title("Energy")
    axes[1].stairs(h_cos[0], h_cos[1])
    axes[1].set_title("cos(theta)")
    axes[2].stairs(h_phi[0], h_phi[1])
    axes[2].set_title("phi")
    axes[3].pcolormesh(x_edges, y_edges, h_hit_map[0].T, cmap="viridis")
    axes[3].set_title("Hit map")
    axes[4].stairs(h_chi2[0], h_chi2[1])
    axes[4].set_title("chi2")
    hit_mean, hit_sigma = _hist_stats(h_resid)
    axes[5].stairs(h_resid[0], h_resid[1])
    axes[5].set_title(f"Hit residuals σ={hit_sigma:.2f} cm")
    axes[6].stairs(h_sipm[0], h_sipm[1])
    axes[6].set_title("SiPM photons")
    axes[7].plot(energy_centers, eff_energy)
    axes[7].set_title("Efficiency vs E")
    axes[8].plot(angle_centers, eff_angle)
    axes[8].set_title("Efficiency vs cos(theta)")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "summary_plot.png", dpi=200)
    plt.close(fig)

    print(f"Plots saved to {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
