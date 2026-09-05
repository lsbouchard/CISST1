#!/usr/bin/env python3
"""Passive contact-resolved spin/charge model; all example inputs are synthetic."""

from dataclasses import dataclass
import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve_banded, solve, eigh
from scipy.optimize import brentq, least_squares

METADATA = {"CreationDate": None, "ModDate": None,
            "Title": "Contact-resolved reciprocal chiral spin transport"}


@dataclass
class Device:
    """Finite-volume weak form C dh/dt = -K h + B (V,hL,hR).

    Physical inputs: length [m], D [m^2/s], gamma [1/s], susceptibility
    [J s^2/m^3], kappa [m/s], a [C/m^3], G [S/m^2]. Examples use
    explicitly dimensionless reference units. I_L and I_R point inward.
    """
    length: float = 1.0
    D: float = 1.0
    gamma: float = 1.0
    susceptibility: float = 1.0
    kappa_left: float = 0.3
    kappa_right: float = 0.3
    a: float = 0.7
    G: float = 1.0
    cells: int = 120
    geometry: str = "lumped"

    def matrices(self):
        if self.geometry not in {"lumped", "local"}:
            raise ValueError("geometry must be lumped or local")
        positive = [self.length, self.D, self.gamma, self.susceptibility, self.G]
        if not np.isfinite(positive + [self.a]).all() or min(positive) <= 0:
            raise ValueError("finite positive length, D, gamma, susceptibility, G required")
        if not isinstance(self.cells, (int, np.integer)) or self.cells < 2:
            raise ValueError("cells must be an integer >= 2")
        contacts = np.asarray([self.kappa_left, self.kappa_right])
        if np.isnan(contacts).any() or (contacts < 0).any():
            raise ValueError("nonnegative contacts required; +inf denotes absorption")
        dx = self.length / self.cells
        c = self.susceptibility * dx
        edge = self.D * self.susceptibility / dx
        # A boundary half-cell and the interface conductance act in series.
        g = np.array([0. if k == 0 else self.susceptibility / (1 / k + dx / (2 * self.D))
                      for k in contacts])
        diag = np.full(self.cells, self.gamma * c + 2 * edge)
        diag[0] += g[0] - edge
        diag[-1] += g[1] - edge
        off = np.full(self.cells - 1, -edge)
        B = np.zeros((self.cells, 3))
        B[:, 0] = self.a * dx
        B[0, 1], B[-1, 2] = g
        return c, diag, off, B, np.array([self.G, *g])

    def stiffness(self):
        """Positive spin-loss matrix, including local charge feedback if selected."""
        c, diag, off, B, direct = self.matrices()
        K = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
        if self.geometry == "local":
            dx = self.length / self.cells
            coefficient = self.a**2 * self.length / self.G
            K += coefficient * (dx*np.eye(self.cells) - dx**2/self.length*np.ones((self.cells, self.cells)))
        return c, K, B, direct

    def response(self, omega):
        if not np.isfinite(omega):
            raise ValueError("omega must be finite")
        c, diag, off, B, direct = self.matrices()
        band = np.zeros((3, self.cells), dtype=complex)
        band[1] = diag - 1j * omega * c
        band[0, 1:] = off
        band[2, :-1] = off
        if self.geometry == "local":
            c, K, B, direct = self.stiffness()
            states = solve(K - 1j*omega*c*np.eye(self.cells), B, assume_a="sym")
        else:
            states = solve_banded((1, 1), band, B)
        E = np.diag([1., -1., -1.])
        Y = np.diag(direct) + E @ B.T @ states
        return Y, states

    def decay_rates(self, load=0.0, count=3):
        if np.isnan(load) or load < 0:
            raise ValueError("nonnegative load required")
        c, K, B, _ = self.stiffness()
        if not isinstance(count, (int, np.integer)) or not 1 <= count <= self.cells:
            raise ValueError("count must be an integer between 1 and cells")
        r = 1 / self.G if np.isinf(load) else load / (1 + load * self.G)
        K += r * np.outer(B[:, 0], B[:, 0])
        return eigh(K / c, subset_by_index=[0, count - 1], eigvals_only=True)


def asymmetric_pole(length, D, gamma, kappa_left, kappa_right):
    """Lowest scalar Robin rate for finite nonnegative, possibly unequal contacts.

    This does not include local electrical feedback. A phase formulation avoids
    multiplying large Biot numbers and covers a single reflecting interface.
    """
    if not np.isfinite([length, D, gamma, kappa_left, kappa_right]).all():
        raise ValueError("finite inputs required; use slowest_pole for two absorbing contacts")
    if min(length, D, gamma) <= 0 or min(kappa_left, kappa_right) < 0:
        raise ValueError("positive length, D, gamma and nonnegative contacts required")
    if kappa_left == kappa_right:
        return slowest_pole(length, D, gamma, kappa_left)
    bl, br = kappa_left*length/D, kappa_right*length/D
    # q = atan(bL/q) + atan(bR/q) is monotone across the fundamental branch.
    q = brentq(lambda q: q - np.arctan2(bl, q) - np.arctan2(br, q),
               0., np.pi, xtol=1e-14)
    return gamma + D*(q/length)**2


def slowest_pole(length, D, gamma, kappa):
    """Exact even Robin root; no drift, identical contacts, positive D/gamma."""
    if not np.isfinite([length, D, gamma]).all() or min(length, D, gamma) <= 0:
        raise ValueError("positive finite length, D, gamma required")
    if np.isnan(kappa) or kappa < 0:
        raise ValueError("nonnegative kappa required")
    if kappa == 0:
        return gamma
    if np.isinf(kappa):
        return gamma + D * (np.pi / length)**2
    b = kappa * length / (2 * D)
    if b <= 1:
        x = brentq(lambda x: x * np.sin(x) - b * np.cos(x), 0., np.pi / 2, xtol=1e-14)
    else:
        # Complementary angle avoids tan overflow near an absorbing contact.
        y = brentq(lambda y: (np.pi / 2 - y) * np.cos(y) / b - np.sin(y),
                   0., np.pi / 2, xtol=1e-14)
        x = np.pi / 2 - y
    return gamma + D * (2 * x / length)**2


def infer_parameters(lengths, rates, relative_sigma=0.02):
    """Fit log(D, gamma, kappa), scaling internally to be independent of SI units."""
    lengths, rates = np.asarray(lengths), np.asarray(rates)
    if lengths.shape != rates.shape or lengths.ndim != 1 or len(lengths) < 3:
        raise ValueError("matching one-dimensional arrays with at least three observations required")
    if not np.isfinite(np.r_[lengths, rates, relative_sigma]).all() or min(lengths.min(), rates.min(), relative_sigma) <= 0:
        raise ValueError("positive finite lengths, rates, and uncertainty required")
    length_scale = np.exp(np.mean(np.log(lengths)))
    rate_scale = np.exp(np.mean(np.log(rates)))
    scaled_lengths, scaled_rates = lengths / length_scale, rates / rate_scale
    def residual(log_parameters):
        D, gamma, kappa = np.exp(log_parameters)
        prediction = [slowest_pole(L, D, gamma, kappa) for L in scaled_lengths]
        return (np.log(prediction) - np.log(scaled_rates)) / relative_sigma
    fit = least_squares(residual, np.log([1., .2, 1.]), bounds=(-12., 12.),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12)
    singular = np.linalg.svd(fit.jac, compute_uv=False)
    if not fit.success or fit.active_mask.any() or singular[-1] <= 1e-8 * singular[0]:
        raise ValueError("fit failed or parameters are locally unidentifiable; broaden the length range")
    covariance = np.linalg.inv(fit.jac.T @ fit.jac)
    fit.parameter_scales = np.array([rate_scale * length_scale**2, rate_scale, rate_scale * length_scale])
    fit.x += np.log(fit.parameter_scales)
    return fit, singular, covariance


def node_rate(device, load):
    """Uniform-density approximation with the physical interface conductances."""
    capacity = device.susceptibility * device.length
    escape = device.susceptibility * (device.kappa_left + device.kappa_right)
    vertex = device.a * device.length
    feedback = 1 / device.G if np.isinf(load) else load / (1 + load * device.G)
    return device.gamma + escape / capacity + vertex**2 * feedback / capacity


def write_csv(name, rows):
    with (ROOT / "data" / name).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    for folder in ["data", "figures", "tables"]:
        (ROOT / folder).mkdir(exist_ok=True)
    device = Device()
    E = np.diag([1., -1., -1.])
    omega = np.logspace(-3, 3, 160)
    rows = []
    for w in omega:
        Y, _ = device.response(w)
        rows.append(dict(omega=w, Y_cL_real=Y[0, 1].real, Y_cL_imag=Y[0, 1].imag,
                         Y_Lc_real=Y[1, 0].real, Y_Lc_imag=Y[1, 0].imag,
                         common_inverse_abs=abs(Y[0, 1] + Y[0, 2]),
                         differential_inverse_abs=abs(Y[0, 1] - Y[0, 2]),
                         reciprocity_error=float(np.max(np.abs(Y - E @ Y.T @ E))),
                         minimum_dissipation_eigenvalue=float(np.linalg.eigvalsh((Y + Y.conj().T) / 2).min())))
    write_csv("device_admittance.csv", rows)
    loads = np.logspace(-3, 3, 90)
    short = device.decay_rates()[0]
    load_rates = np.array([device.decay_rates(R)[0] for R in loads])
    write_csv("device_load_decay.csv", [dict(RG=R*device.G, spatial_rate=rate, short_rate=short,
                                            node_rate=node_rate(device, R))
                                         for R, rate in zip(loads, load_rates)])
    truth = np.array([1., .2, 1.])
    lengths = np.geomspace(.05, 30., 18)
    exact = np.array([slowest_pole(L, *truth) for L in lengths])
    sigma = .02
    observed = exact * np.exp(np.random.default_rng(20260904).normal(0., sigma, len(lengths)))
    fit, singular, covariance = infer_parameters(lengths, observed, sigma)
    parameters = np.exp(fit.x)
    fitted = np.array([slowest_pole(L, *parameters) for L in lengths])
    write_csv("synthetic_thickness.csv", [dict(length=L, exact_rate=v, synthetic_rate=y,
                                              fitted_rate=f, log_standard_deviation=sigma)
                                           for L, v, y, f in zip(lengths, exact, observed, fitted)])
    param_rows = [dict(parameter=name, truth=t, fitted=f,
                       local_log_standard_error=e)
                  for name, t, f, e in zip(["D", "Gamma", "kappa"], truth, parameters, np.sqrt(np.diag(covariance)))]
    write_csv("synthetic_parameter_recovery.csv", param_rows)
    summary = dict(seed=20260904, units="dimensionless synthetic reference units",
                   device_parameters=device.__dict__, truth=dict(zip(["D", "Gamma", "kappa"], truth)),
                   fit_success=bool(fit.success), singular_values=singular.tolist(),
                   fit_parameter_scales=fit.parameter_scales.tolist(),
                   log_covariance=covariance.tolist(), reduced_chi_squared=float(np.sum(fit.fun**2) / (len(lengths) - 3)),
                   short_rate=float(short), open_rate=float(device.decay_rates(np.inf)[0]))
    (ROOT / "data/device_parameters.json").write_text(json.dumps(summary, indent=2) + "\n")
    table = [r"\begin{tabular}{lrrr}", r"\toprule",
             r"Parameter & Input & Inferred & Local log s.e. \\", r"\midrule"]
    for row in param_rows:
        label = {"D": "$D_s$", "Gamma": "$T_1^{-1}$", "kappa": "$\\kappa$"}[row["parameter"]]
        table.append(f"{label} & {row['truth']:.3f} & {row['fitted']:.3f} & {row['local_log_standard_error']:.3f} " + r"\\")
    table += [r"\bottomrule", r"\end{tabular}"]
    (ROOT / "tables/synthetic_inference.tex").write_text("\n".join(table) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
