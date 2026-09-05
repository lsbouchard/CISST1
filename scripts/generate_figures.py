#!/usr/bin/env python3
"""Generate the manuscript figures, source CSV files, and Barnett table."""

from pathlib import Path
import csv
import os

ROOT = Path(__file__).resolve().parents[1]
os.environ["MPLBACKEND"] = "Agg"
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
(ROOT / ".mplconfig").mkdir(exist_ok=True)

import numpy as np
import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import LogFormatterMathtext


FIG = ROOT / "figures"
DATA = ROOT / "data"
TAB = ROOT / "tables"

HBAR = 1.054571817e-34
KB = 1.380649e-23
GAMMA_E_ABS = 1.76085963023e11
N_TURN = 10.5

PDF_METADATA = {
    "Title": "CISST1 generated figure",
    "Author": "Mohamad Niknam and Louis-S. Bouchard",
    "Creator": "scripts/generate_figures.py",
    "CreationDate": None,
    "ModDate": None,
}


def configure_plotting():
    FIG.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    TAB.mkdir(exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "legend.fontsize": 7,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_pdf(fig, filename, *, tight=False):
    kwargs = {"metadata": PDF_METADATA, "dpi": 300}
    if tight:
        kwargs["bbox_inches"] = "tight"
    fig.savefig(FIG / filename, **kwargs)
    plt.close(fig)


def write_csv(filename, header, rows):
    with (DATA / filename).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def sci_tex(value):
    if value == 0:
        return r"$0$"
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / (10**exponent)
    if -2 <= exponent <= 2:
        return f"${value:.3g}$"
    return rf"${mantissa:.3g}\times10^{{{exponent}}}$"


def generate_model_schematic():
    fig, ax = plt.subplots(figsize=(7.0, 2.05))
    ax.set_axis_off()
    boxes = [
        (0.02, 0.37, 0.18, 0.31, "charge current\nthrough chiral axis"),
        (0.26, 0.37, 0.19, 0.31, "current-induced\nelectronic response"),
        (0.50, 0.37, 0.20, 0.31, "SOC-enabled\nspin cross-response"),
        (0.75, 0.34, 0.22, 0.37, "reservoir exchange,\nrelaxation, transport,\nand reciprocal ports"),
    ]
    for x0, y0, width, height, text in boxes:
        patch = FancyBboxPatch(
            (x0, y0),
            width,
            height,
            boxstyle="round,pad=0.012",
            linewidth=1.0,
            facecolor="white",
            edgecolor="black",
        )
        ax.add_patch(patch)
        ax.text(x0 + width / 2, y0 + height / 2, text, ha="center", va="center")

    for left, right in zip(boxes[:-1], boxes[1:]):
        x1 = left[0] + left[2]
        y1 = left[1] + left[3] / 2
        x2 = right[0]
        y2 = right[1] + right[3] / 2
        ax.add_patch(
            FancyArrowPatch(
                (x1 + 0.012, y1),
                (x2 - 0.012, y2),
                arrowstyle="->",
                mutation_scale=12,
                linewidth=1.0,
            )
        )

    ax.text(
        0.11,
        0.15,
        r"$\mathbf{j}_c=j_c\hat{\mathbf{h}}$",
        ha="center",
    )
    ax.text(
        0.355,
        0.15,
        r"$\delta\mathcal{J}^{\rm ang}_h=C_{{\rm ang},j}j_c$",
        ha="center",
        fontsize=7,
    )
    ax.text(
        0.595,
        0.15,
        r"$s_{\chi,h}=\mathcal{A}_{hh}\Omega_{\chi,h}$",
        ha="center",
    )
    ax.text(
        0.86,
        0.15,
        r"$1+D_sq^2T_1+iqv_dT_1-i\omega T_1$",
        ha="center",
        fontsize=6.5,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    save_pdf(fig, "fig_model_schematic.pdf", tight=True)


def generate_barnett_outputs():
    write_csv(
        "model_constants.csv",
        ["constant", "value", "unit"],
        [
            ("hbar", HBAR, "J s"),
            ("kB", KB, "J K^-1"),
            ("abs_gamma_e", GAMMA_E_ABS, "s^-1 T^-1"),
            ("N_turn", N_TURN, "dimensionless"),
        ],
    )

    rates = np.logspace(8, 14, 601)
    omega = 2 * np.pi * rates / N_TURN
    field = omega / GAMMA_E_ABS
    temperatures = (300.0, 100.0, 10.0)
    polarization = {
        temperature: np.tanh(HBAR * omega / (2 * KB * temperature))
        for temperature in temperatures
    }

    write_csv(
        "barnett_estimate.csv",
        ["R_hop_s_inv", "Omega_rad_s_inv", "B_equiv_T"]
        + [f"P_eq_T_{int(temperature)}" for temperature in temperatures],
        (
            [rates[index], omega[index], field[index]]
            + [polarization[temperature][index] for temperature in temperatures]
            for index in range(len(rates))
        ),
    )

    selected_rates = np.array([1e8, 1e9, 1e10, 1e11, 1e12, 1e13])
    selected_omega = 2 * np.pi * selected_rates / N_TURN
    selected_field = selected_omega / GAMMA_E_ABS
    selected_polarization = np.tanh(
        HBAR * selected_omega / (2 * KB * 300.0)
    )
    write_csv(
        "barnett_selected_values.csv",
        ["R_hop_s_inv", "Omega_rad_s_inv", "B_equiv_T", "P_eq_percent_300K"],
        zip(
            selected_rates,
            selected_omega,
            selected_field,
            100 * selected_polarization,
        ),
    )

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.45))
    for temperature in temperatures:
        axes[0].semilogx(
            rates,
            100 * polarization[temperature],
            label=f"{int(temperature)} K",
        )
    marked_rate = 1e13
    marked_value = 100 * np.tanh(
        HBAR * (2 * np.pi * marked_rate / N_TURN) / (2 * KB * 300.0)
    )
    axes[0].plot(marked_rate, marked_value, "o", color="black", markersize=3)
    axes[0].annotate(
        "7.6% at 300 K",
        xy=(marked_rate, marked_value),
        xytext=(2e10, 23),
        arrowprops={"arrowstyle": "->", "linewidth": 0.7},
        fontsize=7,
    )
    axes[0].set_xlabel(r"hopping rate $R_{\rm hop}$ (s$^{-1}$)")
    axes[0].set_ylabel(r"equilibrium $P_{\rm B}$ (%)")
    axes[0].set_ylim(0, 105)
    axes[0].legend(frameon=False, loc="upper left")
    axes[0].set_title("(a) Hypothetical thermal alignment", fontsize=8)

    axes[1].loglog(rates, field)
    axes[1].set_xlabel(r"hopping rate $R_{\rm hop}$ (s$^{-1}$)")
    axes[1].set_ylabel(r"$|B_{\rm B}|$ (T)")
    axes[1].grid(True, which="major", linewidth=0.3, alpha=0.4)
    axes[1].set_title("(b) Trajectory-substituted field", fontsize=8)
    fig.tight_layout(pad=0.5)
    save_pdf(fig, "fig_barnett_estimate.pdf")

    with (TAB / "barnett_table.tex").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("% Generated by scripts/generate_figures.py.\n")
        handle.write(r"\begin{ruledtabular}" + "\n")
        handle.write(r"\begin{tabular}{cccc}" + "\n")
        handle.write(
            r"$R_{\rm hop}$ (s$^{-1}$) & $\Omega_{\rm traj}$ (s$^{-1}$) & "
            r"$|B_{\rm B}|$ (T) & $P_{\rm B}$ at 300 K (\%)\\" + "\n"
        )
        for rate, angular_rate, magnetic_field, spin_polarization in zip(
            selected_rates,
            selected_omega,
            selected_field,
            selected_polarization,
        ):
            handle.write(
                f"{sci_tex(rate)} & {sci_tex(angular_rate)} & "
                f"{sci_tex(magnetic_field)} & "
                f"{sci_tex(100 * spin_polarization)}"
                + r"\\"
                + "\n"
            )
        handle.write(r"\end{tabular}" + "\n")
        handle.write(r"\end{ruledtabular}" + "\n")


def generate_relaxation_response():
    omega_t1 = np.logspace(-2, 2, 800)
    amplitude = 1 / np.sqrt(1 + omega_t1**2)
    phase_deg = np.degrees(np.arctan(omega_t1))
    write_csv(
        "relaxation_response.csv",
        ["omega_T1", "normalized_amplitude", "lag_deg"],
        zip(omega_t1, amplitude, phase_deg),
    )

    fig, left_axis = plt.subplots(figsize=(3.5, 2.45))
    left_axis.semilogx(omega_t1, amplitude, label="amplitude")
    left_axis.set_xlabel(r"normalized frequency $\omega T_1$")
    left_axis.set_ylabel(r"normalized amplitude")
    left_axis.set_ylim(0, 1.05)
    left_axis.grid(True, which="major", linewidth=0.3, alpha=0.4)
    right_axis = left_axis.twinx()
    right_axis.semilogx(omega_t1, phase_deg, linestyle="--", label="lag")
    right_axis.set_ylabel(r"phase lag (degrees)")
    right_axis.set_ylim(0, 90)
    lines_left, labels_left = left_axis.get_legend_handles_labels()
    lines_right, labels_right = right_axis.get_legend_handles_labels()
    left_axis.legend(
        lines_left + lines_right,
        labels_left + labels_right,
        frameon=False,
        loc="center left",
    )
    fig.tight_layout(pad=0.5)
    save_pdf(fig, "fig_relaxation_response.pdf")


def generate_length_scaling():
    normalized_length = np.linspace(0, 6, 700)
    drift_density = 1 - np.exp(-normalized_length)
    diffusive_density = 1 - 1 / np.cosh(normalized_length)
    outward_flux = np.tanh(normalized_length)
    write_csv(
        "length_scaling.csv",
        [
            "normalized_length",
            "drift_exit_density_over_s_star",
            "diffusive_exit_density_over_s_star",
            "left_outward_flux_over_Ds_s_star_per_ell_s",
        ],
        zip(
            normalized_length,
            drift_density,
            diffusive_density,
            outward_flux,
        ),
    )

    fig, ax = plt.subplots(figsize=(3.5, 2.45))
    ax.plot(
        normalized_length,
        drift_density,
        label=r"drift density: $1-e^{-L/\ell_d}$",
    )
    ax.plot(
        normalized_length,
        diffusive_density,
        linestyle="--",
        label=r"diffusive density: $1-{\rm sech}(L/\ell_s)$",
    )
    ax.plot(
        normalized_length,
        outward_flux,
        linestyle=":",
        label=r"outward flux: $\tanh(L/\ell_s)$",
    )
    ax.set_xlabel(r"dimensionless thickness ($L/\ell_d$ or $L/\ell_s$)")
    ax.set_ylabel(r"normalized response")
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 1.05)
    ax.grid(True, linewidth=0.3)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout(pad=0.5)
    save_pdf(fig, "fig_length_scaling.pdf")


def generate_diffusion_landscape():
    t1_ns = np.logspace(-2, 2, 200)
    diffusion = np.logspace(-9, -5, 220)
    t1_grid, diffusion_grid = np.meshgrid(t1_ns, diffusion)
    ell_nm = np.sqrt(diffusion_grid * t1_grid * 1e-9) * 1e9

    write_csv(
        "diffusion_landscape_grid.csv",
        ["T1_ns", "D_s_m2_s", "ell_s_nm"],
        (
            [t1_ns[column], diffusion[row], ell_nm[row, column]]
            for row in range(len(diffusion))
            for column in range(len(t1_ns))
        ),
    )

    fig, ax = plt.subplots(figsize=(3.7, 2.7))
    color_mesh = ax.pcolormesh(
        t1_ns,
        diffusion,
        ell_nm,
        shading="auto",
        norm=LogNorm(vmin=0.1, vmax=1000),
        cmap="viridis",
        edgecolors="face",
        antialiased=False,
        rasterized=True,
        snap=True,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$T_1$ (ns)")
    ax.set_ylabel(r"$D_s$ (m$^2$ s$^{-1}$)")
    ax.set_title("Kinematic map; continuum validity is not implied", fontsize=8)
    colorbar = fig.colorbar(
        color_mesh,
        ax=ax,
        ticks=[0.1, 1, 10, 100, 1000],
        format=LogFormatterMathtext(),
    )
    colorbar.set_label(r"$\ell_s=\sqrt{D_sT_1}$ (nm)")
    contours = ax.contour(
        t1_grid,
        diffusion_grid,
        ell_nm,
        levels=[0.3, 1, 3, 10, 30, 100, 300],
        colors="black",
        linewidths=0.4,
    )
    ax.clabel(contours, fmt="%g", fontsize=6)
    fig.tight_layout(pad=0.5)
    save_pdf(fig, "fig_spin_diffusion_landscape.pdf")


def main():
    configure_plotting()
    generate_model_schematic()
    generate_barnett_outputs()
    generate_relaxation_response()
    generate_length_scaling()
    generate_diffusion_landscape()
    print("Wrote figures, data, and table to", ROOT)


if __name__ == "__main__":
    main()
