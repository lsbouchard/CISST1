#!/usr/bin/env python3
"""Analytic finite-contact solution for the 1D relaxation-diffusion model."""

from dataclasses import dataclass
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


FIG = ROOT / "figures"
DATA = ROOT / "data"
PDF_METADATA = {
    "Title": "Finite contact transparency",
    "Author": "Mohamad Niknam and Louis-S. Bouchard",
    "Creator": "scripts/solve_contact_transparency.py",
    "CreationDate": None,
    "ModDate": None,
}


@dataclass(frozen=True)
class AnalyticProfile:
    """Solution u=1+A exp(r_plus x)+B exp(r_minus x)."""

    length: float
    peclet: float
    r_plus: float
    r_minus: float
    coefficient_plus: float
    coefficient_minus: float

    def value(self, x):
        x_array = np.asarray(x)
        return (
            1.0
            + self.coefficient_plus * np.exp(self.r_plus * x_array)
            + self.coefficient_minus * np.exp(self.r_minus * x_array)
        )

    def derivative(self, x):
        x_array = np.asarray(x)
        return (
            self.coefficient_plus
            * self.r_plus
            * np.exp(self.r_plus * x_array)
            + self.coefficient_minus
            * self.r_minus
            * np.exp(self.r_minus * x_array)
        )

    def second_derivative(self, x):
        x_array = np.asarray(x)
        return (
            self.coefficient_plus
            * self.r_plus**2
            * np.exp(self.r_plus * x_array)
            + self.coefficient_minus
            * self.r_minus**2
            * np.exp(self.r_minus * x_array)
        )

    def average(self):
        plus_integral = np.expm1(self.r_plus * self.length) / self.r_plus
        minus_integral = np.expm1(self.r_minus * self.length) / self.r_minus
        return (
            self.length
            + self.coefficient_plus * plus_integral
            + self.coefficient_minus * minus_integral
        ) / self.length


def solve_profile(length_over_ell, tau, peclet=0.0, absorbing=False):
    """Solve u''-Pe*u'-(u-1)=0 with symmetric outward-flux contacts.

    The dimensionless transparency is tau=kappa*ell_s/D_s. Setting
    absorbing=True imposes u=0 at both contacts, the tau -> infinity limit.
    """

    if length_over_ell <= 0:
        raise ValueError("length_over_ell must be positive")
    if not absorbing and tau < 0:
        raise ValueError("tau must be nonnegative")

    discriminant = np.sqrt(peclet**2 + 4.0)
    r_plus = 0.5 * (peclet + discriminant)
    r_minus = 0.5 * (peclet - discriminant)
    length = float(length_over_ell)
    exp_plus = np.exp(r_plus * length)
    exp_minus = np.exp(r_minus * length)

    if absorbing:
        matrix = np.array(
            [
                [1.0, 1.0],
                [exp_plus, exp_minus],
            ]
        )
        rhs = np.array([-1.0, -1.0])
    else:
        matrix = np.array(
            [
                [
                    r_plus - peclet - tau,
                    r_minus - peclet - tau,
                ],
                [
                    exp_plus * (r_plus - peclet + tau),
                    exp_minus * (r_minus - peclet + tau),
                ],
            ]
        )
        rhs = np.array([peclet + tau, peclet - tau])

    coefficient_plus, coefficient_minus = np.linalg.solve(matrix, rhs)
    return AnalyticProfile(
        length=length,
        peclet=float(peclet),
        r_plus=float(r_plus),
        r_minus=float(r_minus),
        coefficient_plus=float(coefficient_plus),
        coefficient_minus=float(coefficient_minus),
    )


def contact_observables(profile):
    """Return mean density and exact dimensionless outward boundary fluxes."""

    u_left = float(profile.value(0.0))
    u_right = float(profile.value(profile.length))
    du_left = float(profile.derivative(0.0))
    du_right = float(profile.derivative(profile.length))
    left_flux = du_left - profile.peclet * u_left
    right_flux = profile.peclet * u_right - du_right
    return profile.average(), u_left, u_right, left_flux, right_flux


def residuals(profile, tau, absorbing):
    """Return maximum ODE and boundary residuals for an analytic profile."""

    x = np.linspace(0.0, profile.length, 17)
    u = profile.value(x)
    du = profile.derivative(x)
    d2u = profile.second_derivative(x)
    ode_residual = d2u - profile.peclet * du - (u - 1.0)

    if absorbing:
        boundary_residuals = np.array(
            [profile.value(0.0), profile.value(profile.length)]
        )
    else:
        boundary_residuals = np.array(
            [
                profile.derivative(0.0)
                - (profile.peclet + tau) * profile.value(0.0),
                profile.derivative(profile.length)
                - (profile.peclet - tau)
                * profile.value(profile.length),
            ]
        )
    return (
        float(np.max(np.abs(ode_residual))),
        float(np.max(np.abs(boundary_residuals))),
    )


def main():
    FIG.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    peclet = 0.0
    lengths = np.linspace(0.05, 8.0, 220)
    cases = [
        ("absorbing", np.inf, True, "absorbing contacts"),
        ("partial", 1.0, False, r"partly transparent, $\kappa\ell_s/D_s=1$"),
        ("reflecting", 0.0, False, "reflecting contacts"),
    ]

    rows = []
    validation_rows = []
    for case_key, tau, absorbing, label in cases:
        case_ode_residual = 0.0
        case_boundary_residual = 0.0
        for length in lengths:
            tau_solve = 0.0 if absorbing else tau
            profile = solve_profile(
                length,
                tau_solve,
                peclet=peclet,
                absorbing=absorbing,
            )
            average, u_left, u_right, left_flux, right_flux = (
                contact_observables(profile)
            )
            ode_residual, boundary_residual = residuals(
                profile,
                tau_solve,
                absorbing,
            )
            case_ode_residual = max(case_ode_residual, ode_residual)
            case_boundary_residual = max(
                case_boundary_residual,
                boundary_residual,
            )
            rows.append(
                {
                    "case": case_key,
                    "label": label,
                    "L_over_ell_s": length,
                    "tau_left": "inf" if absorbing else tau,
                    "tau_right": "inf" if absorbing else tau,
                    "peclet": peclet,
                    "average_s_over_s_star": average,
                    "left_boundary_s_over_s_star": u_left,
                    "right_boundary_s_over_s_star": u_right,
                    "left_outward_flux_norm": left_flux,
                    "right_outward_flux_norm": right_flux,
                    "total_outward_flux_norm": left_flux + right_flux,
                }
            )
        validation_rows.append(
            {
                "case": case_key,
                "max_abs_ode_residual": case_ode_residual,
                "max_abs_boundary_residual": case_boundary_residual,
            }
        )

    for validation in validation_rows:
        if validation["max_abs_ode_residual"] > 1e-10:
            raise RuntimeError(f"ODE validation failed: {validation}")
        if validation["max_abs_boundary_residual"] > 1e-10:
            raise RuntimeError(f"boundary validation failed: {validation}")

    output_csv = DATA / "contact_transparency.csv"
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    validation_csv = DATA / "contact_validation.csv"
    with validation_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(validation_rows[0].keys()),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(validation_rows)

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

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.55))
    colors = {
        "absorbing": "#1f77b4",
        "partial": "#d62728",
        "reflecting": "#2ca02c",
    }
    styles = {"absorbing": "-", "partial": "--", "reflecting": ":"}
    for case_key, _tau, _absorbing, label in cases:
        case_rows = [row for row in rows if row["case"] == case_key]
        x_values = np.array([row["L_over_ell_s"] for row in case_rows])
        average = np.array(
            [row["average_s_over_s_star"] for row in case_rows]
        )
        flux_per_contact = (
            np.array(
                [row["total_outward_flux_norm"] for row in case_rows]
            )
            / 2.0
        )
        axes[0].plot(
            x_values,
            average,
            styles[case_key],
            color=colors[case_key],
            label=label,
        )
        axes[1].plot(
            x_values,
            flux_per_contact,
            styles[case_key],
            color=colors[case_key],
            label=label,
        )

    axes[0].set_xlabel(r"layer thickness $L/\ell_s$")
    axes[0].set_ylabel(r"mean density $\langle s\rangle/s_*$")
    axes[0].set_xlim(0, 8)
    axes[0].set_ylim(0, 1.05)
    axes[0].grid(True, linewidth=0.3)
    axes[0].text(
        0.03,
        0.95,
        "(a) accumulation",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
    )

    axes[1].set_xlabel(r"layer thickness $L/\ell_s$")
    axes[1].set_ylabel(
        r"outward flux per contact $\Phi/(D_ss_*/\ell_s)$"
    )
    axes[1].set_xlim(0, 8)
    axes[1].set_ylim(0, 1.05)
    axes[1].grid(True, linewidth=0.3)
    axes[1].text(
        0.03,
        0.95,
        "(b) contact loss",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
    )
    axes[1].legend(frameon=False, loc="lower right")

    fig.tight_layout(pad=0.5)
    fig.savefig(
        FIG / "fig_contact_transparency.pdf",
        metadata=PDF_METADATA,
    )
    plt.close(fig)

    print("Wrote", output_csv, validation_csv, "and contact figure")


if __name__ == "__main__":
    main()
