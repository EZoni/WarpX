#!/usr/bin/env python3

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as scc
from openpmd_viewer import OpenPMDTimeSeries

REACTION_CONFIG = {
    "DD": {
        "e_min": 0.0,
        "e_max": 15.0,
        "output": "dd_anisotropic_neutron_spectrum.png",
        "test_name": "test_3d_deuterium_deuterium_fusion_anisotropic_beam_target",
        "xticks": [0.0, 5.0, 10.0, 15.0],
        "scale": 1.0,
    },
    "DT": {
        "e_min": 10.0,
        "e_max": 30.0,
        "output": "dt_anisotropic_neutron_spectrum.png",
        "test_name": "test_3d_deuterium_tritium_fusion_anisotropic_beam_target",
        "xticks": [10.0, 15.0, 20.0, 25.0, 30.0],
        "scale": 2.6,
    },
}
BEAM_SPEEDS_PERCENT = [1, 5, 10]


def parse_reaction(value):
    reaction = value.upper()
    if reaction not in REACTION_CONFIG:
        choices = ", ".join(REACTION_CONFIG)
        raise argparse.ArgumentTypeError(f"reaction must be one of: {choices}")
    return reaction


def parse_args():
    argv = sys.argv[1:]
    if len(argv) == 2 and argv[1].lower().endswith((".png", ".pdf", ".svg")):
        argv = [argv[0], "--output", argv[1]]

    parser = argparse.ArgumentParser()
    parser.add_argument("diag_dirs", nargs="*")
    parser.add_argument("--reaction", type=parse_reaction, default="DD")
    parser.add_argument("-o", "--output")
    parser.add_argument("--labels", nargs="+")
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = REACTION_CONFIG[args.reaction]["output"]

    return args


def discover_diag_dirs(config):
    """Find completed outputs from this reaction's sibling CTest runs."""
    run_root = Path.cwd().parent
    test_name = config["test_name"]
    diag_dirs = []
    for speed in BEAM_SPEEDS_PERCENT:
        diag_dir = run_root / f"{test_name}_{speed}pc" / "diags" / "diag1"
        if diag_dir.is_dir() and any(diag_dir.glob("openpmd_*")):
            diag_dirs.append(diag_dir)

    # Preserve the command-line script's original single-run behavior outside
    # a CTest run directory.
    if not diag_dirs:
        diag_dirs = [Path("diags/diag1")]
    return diag_dirs


def infer_label(diag_dir):
    path = Path(diag_dir)
    # Prefer the beam velocity encoded in directory names such as run_5pc.
    for part in reversed(path.parts):
        match = re.search(r"(\d+)pc", part)
        if match:
            return rf"$u/c = {match.group(1)}\%$"

    # Fall back to the input deck recorded with the diagnostics.
    for input_path in [
        path / "warpx_used_inputs",
        path.parent / "warpx_used_inputs",
        path.parent.parent / "warpx_used_inputs",
    ]:
        if not input_path.exists():
            continue
        text = input_path.read_text()
        match = re.search(
            r"deuterium_beam\.momentum_function_uz\(x,y,z\)\s*=\s*"
            r"([+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)",
            text,
        )
        if match:
            percent = 100.0 * float(match.group(1))
            return rf"$u/c = {percent:g}\%$"

    return path.parent.name if path.name == "diag1" else path.name


def get_neutron_spectrum(diag_dir):
    ts = OpenPMDTimeSeries(diag_dir)
    iteration = ts.iterations[-1]
    ux, uy, uz, w = ts.get_particle(
        ["ux", "uy", "uz", "w"], species="neutron_1", iteration=iteration
    )

    u2 = ux**2 + uy**2 + uz**2
    u = np.sqrt(u2)
    gamma = np.sqrt(1.0 + u2)
    m_neutron = scc.m_n
    # Convert neutron momentum u = gamma beta to kinetic energy.
    energy_MeV = (gamma - 1.0) * m_neutron * scc.c**2 / scc.e / 1.0e6
    # Measure the emission angle relative to the beam axis (+z).
    cos_theta = np.divide(uz, u, out=np.zeros_like(uz), where=u > 0.0)
    theta_deg = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
    return energy_MeV, theta_deg, w


def plot_neutron_spectra(diag_dirs, labels, config, output):
    e_min = config["e_min"]
    e_max = config["e_max"]
    energy_bins = np.linspace(e_min, e_max, 240)
    # Use the same energy range for the energy-resolved mean emission angle.
    energy_bins_angle = np.linspace(e_min, e_max, 180)
    angle_centers = 0.5 * (energy_bins_angle[:-1] + energy_bins_angle[1:])

    fig, ax_spectrum = plt.subplots(figsize=(6.0, 2.5), constrained_layout=True)
    ax_angle = ax_spectrum.twinx()

    for diag_dir, label in zip(diag_dirs, labels):
        energy_MeV, theta_deg, w = get_neutron_spectrum(diag_dir)
        assert energy_MeV.size > 0, f"No neutron macroparticles found in {diag_dir}."

        # Histogram macroparticle weights to obtain the neutron energy spectrum.
        hist, edges = np.histogram(energy_MeV, bins=energy_bins, weights=w)
        assert hist.sum() > 0.0, (
            f"No neutron weight in the plotted energy range for {diag_dir}."
        )
        bin_width = edges[1] - edges[0]
        # Normalize by total weight and bin width so the plotted area is unity.
        spectrum = hist / (hist.sum() * bin_width)
        # Apply the reaction-dependent ad-hoc scaling to the spectrum.
        spectrum *= config.get("scale", 1.0)
        # Plot against bin centers instead of left bin edges.
        centers = 0.5 * (edges[:-1] + edges[1:])

        (line,) = ax_spectrum.plot(centers, spectrum, lw=1.4, label=label)

        # Overlay the weighted mean emission angle in each neutron-energy bin.
        # The first histogram is the numerator, sum(w * theta).
        weighted_angle, _ = np.histogram(
            energy_MeV, bins=energy_bins_angle, weights=w * theta_deg
        )
        # The second histogram is the denominator, sum(w).
        angle_weight, _ = np.histogram(energy_MeV, bins=energy_bins_angle, weights=w)
        valid = angle_weight > 0.0
        mean_angle = np.empty_like(weighted_angle)
        mean_angle[:] = np.nan
        # Leave empty bins as NaN so Matplotlib breaks the dashed curve there.
        mean_angle[valid] = weighted_angle[valid] / angle_weight[valid]
        ax_angle.plot(
            angle_centers,
            mean_angle,
            color=line.get_color(),
            ls="--",
            lw=1.2,
        )

    ax_spectrum.set_xlim(e_min, e_max)
    ax_spectrum.set_ylim(0.0, 2.2)
    ax_angle.set_ylim(0.0, 180.0)
    ax_spectrum.set_xticks(config["xticks"])
    ax_spectrum.set_yticks([0.0, 1.0, 2.0])
    ax_angle.set_yticks([0.0, 60.0, 120.0, 180.0])
    ax_spectrum.set_xlabel(r"$E_n$ (MeV)")
    ax_spectrum.set_ylabel(r"$dN/dE_n$ (arb. units)")
    ax_angle.set_ylabel(r"$\theta$ (degrees)")
    ax_spectrum.legend(frameon=False, loc="upper right", fontsize=8)
    ax_spectrum.minorticks_on()
    ax_angle.minorticks_on()

    fig.savefig(output, dpi=200)
    print(f"Saved {output}")


def main():
    args = parse_args()
    if not args.diag_dirs:
        args.diag_dirs = discover_diag_dirs(REACTION_CONFIG[args.reaction])

    labels = args.labels
    if labels is None:
        labels = [infer_label(diag_dir) for diag_dir in args.diag_dirs]
    assert len(labels) == len(args.diag_dirs), (
        "Number of labels must match diagnostics."
    )

    plot_neutron_spectra(
        args.diag_dirs,
        labels,
        REACTION_CONFIG[args.reaction],
        args.output,
    )


if __name__ == "__main__":
    main()
