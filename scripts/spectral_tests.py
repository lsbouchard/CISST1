#!/usr/bin/env python3
"""Electrical spectral weights and certified fast-mode truncation bounds.

Inputs use the SI conventions of reciprocal_device.Device. Figure inputs are
synthetic reference units. Bounds exclude measurement and model uncertainty.
"""

import json
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

from reciprocal_device import Device, ROOT, METADATA, plt, write_csv


def spectrum(device):
    """Return rates [1/s] and modal port overlaps (N-by-3, mixed SI units)."""
    capacity, loss, ports, _ = device.stiffness()
    rates, vectors = eigh(loss / capacity)
    overlaps = vectors.T @ ports / np.sqrt(capacity)
    return rates, overlaps


def tail_bracket(rate, weight, dc_tail, cutoff, alpha):
    """Bracket first visible loaded rate from one mode and a gapped dc tail.

    rate, cutoff [1/s]; weight [S/m^2/s]; dc_tail [S/m^2]; alpha [ohm m^2].
    The caller must establish that all omitted visible rates exceed cutoff.
    Raises ValueError when the sufficient gap condition cannot be certified.
    """
    values = np.asarray([rate, weight, dc_tail, cutoff, alpha])
    if not np.isfinite(values).all() or min(rate, weight) <= 0:
        raise ValueError("finite positive first rate/weight required")
    if min(dc_tail, alpha) < 0 or cutoff <= rate:
        raise ValueError("nonnegative dc tail/load and cutoff above first rate required")
    if alpha == 0:
        return rate, rate
    shift = alpha * weight / (1 + alpha * dc_tail)
    upper = rate + shift
    if upper >= cutoff:
        raise ValueError("unresolved-mode gap does not certify this load")
    if dc_tail == 0:
        return upper, upper
    # Solve for a dimensionless shift fraction to avoid pole singularities and
    # preserve accuracy under changes of physical rate and conductance units.
    def equation(fraction):
        gap = 1 - (rate + fraction * shift) / cutoff
        return fraction * gap - (1 + alpha * dc_tail) * gap + alpha * dc_tail * fraction
    fraction = brentq(equation, 0., 1., xtol=1e-14)
    return rate + fraction * shift, upper


def series_correction(widths, resistivity, beta):
    """Piecewise-constant nonuniform local electrical weak form.

    widths [m], resistivity [ohm m], beta [C/m^2], all length N.
    Returns bare G [S/m^2], bulk charge column [C/m^2], and
    positive-semidefinite local loss matrix [J s/m^2].
    """
    widths, resistivity, beta = [np.asarray(x, dtype=float) for x in (widths, resistivity, beta)]
    if widths.ndim != 1 or len(widths) == 0 or widths.shape != resistivity.shape or widths.shape != beta.shape:
        raise ValueError("matching nonempty one-dimensional material arrays required")
    if not np.isfinite(np.r_[widths, resistivity, beta]).all() or min(widths.min(), resistivity.min()) <= 0:
        raise ValueError("positive finite cell widths/resistivity and finite beta required")
    measure = widths * resistivity
    conductance = 1 / measure.sum()
    vertex = conductance * measure * beta
    loss = np.diag(measure * beta**2) - np.outer(vertex, vertex) / conductance
    return conductance, vertex, loss


def example_data(cells=120):
    device = Device(D=.1, gamma=.2, a=1.2, kappa_left=1., kappa_right=.1, cells=cells)
    rates, overlaps = spectrum(device)
    weights = overlaps[:, 0]**2
    dc_tail = float(np.sum(weights[1:] / rates[1:]))
    cutoff = 2.2
    if rates[1] < cutoff:
        raise AssertionError("synthetic omitted spectrum violates chosen cutoff")
    modes = [dict(mode=n+1, rate=rate, charge_weight=weight,
                  left_overlap=overlap[1], right_overlap=overlap[2])
             for n, (rate, weight, overlap) in enumerate(zip(rates, weights, overlaps))]
    frequency = []
    for omega in np.geomspace(.01, 100., 200):
        excess = np.sum(weights / (rates - 1j*omega))
        first = weights[0] / (rates[0] - 1j*omega)
        frequency.append(dict(omega=omega, excess_real=excess.real, excess_imag=excess.imag,
                              first_mode_real=first.real, first_mode_imag=first.imag))
    loads = []
    for resistance in np.geomspace(.001, 100., 100):
        alpha = resistance / (1 + resistance * device.G)
        lower, upper = tail_bracket(rates[0], weights[0], dc_tail, cutoff, alpha)
        exact = device.decay_rates(resistance)[0]
        if not lower-1e-10 <= exact <= upper+1e-10:
            raise AssertionError("independent loaded diagonalization lies outside tail bounds")
        loads.append(dict(RG=resistance*device.G, exact_rate=exact, lower_rate=lower,
                          upper_rate=upper, omitted_tail_rate=rates[0]+alpha*weights[0]))
    summary = dict(units="synthetic dimensionless reference units", parameters=device.__dict__,
                   first_rate=float(rates[0]), first_weight=float(weights[0]),
                   residual_dc=dc_tail, next_rate=float(rates[1]), cutoff=cutoff,
                   total_dc=float(np.sum(weights/rates)),
                   open_rate=float(device.decay_rates(np.inf)[0]))
    return modes, frequency, loads, summary


def main():
    for folder in ("data", "figures"):
        (ROOT/folder).mkdir(exist_ok=True)
    modes, frequency, loads, summary = example_data()
    write_csv("spectral_modes.csv", modes)
    write_csv("spectral_admittance.csv", frequency)
    write_csv("spectral_load_bounds.csv", loads)
    (ROOT/"data/spectral_parameters.json").write_text(json.dumps(summary, indent=2)+"\n")
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.9), layout="constrained")
    omega = [r["omega"] for r in frequency]
    for part, color, label in [("real", "#245e91", "Real"), ("imag", "#ae3e37", "Imaginary")]:
        axes[0].loglog(omega, [r[f"excess_{part}"] for r in frequency], color=color, label=label)
        axes[0].loglog(omega, [r[f"first_mode_{part}"] for r in frequency], color=color, ls="--", lw=1)
    axes[0].plot([], [], color="0.35", ls="--", label="First mode only")
    axes[0].set(xlabel=r"Frequency $\omega$ (reference units)", ylabel=r"Excess admittance $Y_{cc}-G$",
                title="(a) Electrical calibration")
    axes[0].legend(fontsize=8, frameon=False)
    load = [r["RG"] for r in loads]
    axes[1].fill_between(load, [r["lower_rate"] for r in loads], [r["upper_rate"] for r in loads],
                         color="#a4d9bf", label="Certified fast-mode bounds")
    axes[1].semilogx(load, [r["exact_rate"] for r in loads], color="#202020", label="Full device")
    axes[1].semilogx(load, [r["omitted_tail_rate"] for r in loads], color="#8a4c91", ls=":",
                    label="Fast modes discarded")
    axes[1].set(xlabel=r"Electrical load $RG$", ylabel="Loaded rate (reference units)",
                title="(b) Prediction from one mode")
    axes[1].legend(fontsize=8, frameon=False, loc="upper left")
    fig.savefig(ROOT/"figures/fig_spectral_load.pdf", metadata=METADATA)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
