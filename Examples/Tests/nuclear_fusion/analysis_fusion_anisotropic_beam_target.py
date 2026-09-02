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
        "output": "deuterium_deuterium_fusion_anisotropic_beam_target_neutron_spectrum.png",
        "xticks": [0.0, 5.0, 10.0, 15.0],
        "scale": 1.0,
    },
    "DT": {
        "e_min": 10.0,
        "e_max": 30.0,
        "output": "deuterium_tritium_fusion_anisotropic_beam_target_neutron_spectrum.png",
        "xticks": [10.0, 15.0, 20.0, 25.0, 30.0],
        "scale": 2.6,
    },
}

# Reference values obtained with the fixed random seed in the input files. These
# are integral observables, with intentionally broad tolerances for portability
# across compute backends and particle precisions.
REFERENCE_METRICS = {
    "DD": {
        1: {"mean_energy": 2.4961, "std_energy": 0.21466, "mean_cos_theta": 0.04312},
        5: {"mean_energy": 3.6232, "std_energy": 1.3961, "mean_cos_theta": 0.15298},
        10: {"mean_energy": 7.1341, "std_energy": 3.8885, "mean_cos_theta": 0.21013},
    },
    "DT": {
        1: {"mean_energy": 14.0825, "std_energy": 0.37748, "mean_cos_theta": 0.01812},
        5: {"mean_energy": 15.5449, "std_energy": 1.9861, "mean_cos_theta": 0.13423},
        10: {"mean_energy": 19.9295, "std_energy": 4.8190, "mean_cos_theta": 0.20604},
    },
}


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
    parser.add_argument(
        "--plot",
        action="store_true",
        help="plot all supplied diagnostics after validating them",
    )
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = REACTION_CONFIG[args.reaction]["output"]
    else:
        args.plot = True

    return args


def find_used_inputs(diag_dir):
    path = Path(diag_dir)
    for input_path in [
        path / "warpx_used_inputs",
        path.parent / "warpx_used_inputs",
        path.parent.parent / "warpx_used_inputs",
    ]:
        if input_path.exists():
            return input_path
    raise AssertionError(f"Could not find warpx_used_inputs for {diag_dir}.")


def infer_beam_speed_percent(diag_dir):
    text = find_used_inputs(diag_dir).read_text()
    match = re.search(
        r"deuterium_beam\.momentum_function_uz\(x,y,z\)\s*=\s*"
        r"([+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)",
        text,
    )
    assert match, f"Could not determine the beam momentum for {diag_dir}."
    speed_percent = 100.0 * float(match.group(1))
    rounded_speed = round(speed_percent)
    assert np.isclose(speed_percent, rounded_speed), (
        f"Unsupported beam momentum u/c = {speed_percent:g}% in {diag_dir}."
    )
    return rounded_speed


def infer_label(diag_dir):
    speed_percent = infer_beam_speed_percent(diag_dir)
    return rf"$u/c = {speed_percent:g}\%$"


def weighted_mean(values, weights):
    return np.sum(weights * values) / np.sum(weights)


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


def validate_neutron_spectrum(diag_dir, reaction):
    energy_MeV, theta_deg, w = get_neutron_spectrum(diag_dir)
    assert energy_MeV.size >= 1000, (
        f"Too few neutron macroparticles ({energy_MeV.size}) in {diag_dir}."
    )
    assert np.all(np.isfinite(energy_MeV))
    assert np.all(np.isfinite(theta_deg))
    assert np.all(np.isfinite(w)) and np.all(w > 0.0)

    cos_theta = np.cos(np.radians(theta_deg))
    mean_energy = weighted_mean(energy_MeV, w)
    std_energy = np.sqrt(weighted_mean((energy_MeV - mean_energy) ** 2, w))
    mean_cos_theta = weighted_mean(cos_theta, w)
    metrics = {
        "mean_energy": mean_energy,
        "std_energy": std_energy,
        "mean_cos_theta": mean_cos_theta,
    }

    speed_percent = infer_beam_speed_percent(diag_dir)
    assert speed_percent in REFERENCE_METRICS[reaction], (
        f"No reference metrics for u/c = {speed_percent:g}%."
    )
    reference = REFERENCE_METRICS[reaction][speed_percent]
    np.testing.assert_allclose(
        [metrics["mean_energy"], metrics["std_energy"]],
        [reference["mean_energy"], reference["std_energy"]],
        rtol=0.05,
        atol=0.0,
    )
    np.testing.assert_allclose(
        metrics["mean_cos_theta"],
        reference["mean_cos_theta"],
        rtol=0.0,
        atol=0.03,
    )
    print(
        f"{reaction}, u/c = {speed_percent:g}%: "
        f"mean energy = {mean_energy:.6g} MeV, "
        f"energy std. dev. = {std_energy:.6g} MeV, "
        f"mean cos(theta) = {mean_cos_theta:.6g}"
    )


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
        args.diag_dirs = [Path("diags/diag1")]

    for diag_dir in args.diag_dirs:
        validate_neutron_spectrum(diag_dir, args.reaction)

    if args.plot:
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
