#!/usr/bin/env python3
"""Model-misspecification, dynamic-reduction, and electrical-load calculations.

All parameter sets are synthetic reference units, not experimental estimates.
"""
import json
import numpy as np
from scipy.linalg import eigh
from reciprocal_device import (ROOT, METADATA, Device, asymmetric_pole,
                               slowest_pole, infer_parameters, node_rate, write_csv, load_feedback)
import matplotlib.pyplot as plt


def contact_fits():
    samples, fits = [], []
    for low, high in [(.05, 30.), (.5, 5.), (1., 3.)]:
        lengths = np.geomspace(low, high, 18)
        exact = np.array([asymmetric_pole(L, 1., .2, 1., .1) for L in lengths])
        fit, _, _ = infer_parameters(lengths, exact)
        parameters = np.exp(fit.x)
        predicted = np.array([slowest_pole(L, *parameters) for L in lengths])
        residual = np.log(predicted/exact)
        fits.append(dict(length_min=low, length_max=high, D=parameters[0],
                         Gamma=parameters[1], kappa=parameters[2],
                         rms_log_mismatch=np.sqrt(np.mean(residual**2))))
        samples.extend(dict(length_min=low, length_max=high, length=L,
                            asymmetric_rate=y, fitted_rate=f, log_residual=r)
                       for L, y, f, r in zip(lengths, exact, predicted, residual))
    return samples, fits


def length_information():
    rows = []
    theta = np.array([1., .2, 1.])
    step = 1e-4
    for low, high in [(.05, 30.), (.5, 5.), (1., 3.), (.05, .2), (10., 30.)]:
        lengths = np.geomspace(low, high, 18)
        J = np.empty((18, 3))
        for k in range(3):
            p, m = theta.copy(), theta.copy()
            p[k] *= np.exp(step)
            m[k] *= np.exp(-step)
            J[:, k] = (np.log([slowest_pole(L, *p) for L in lengths])
                       - np.log([slowest_pole(L, *m) for L in lengths]))/(2*step*.02)
        _, singular, vt = np.linalg.svd(J, full_matrices=False)
        covariance = (vt.T/singular**2) @ vt
        errors = np.sqrt(np.diag(covariance))
        rows.append(dict(length_min=low, length_max=high, log_se_D=errors[0],
                         log_se_Gamma=errors[1], log_se_kappa=errors[2],
                         condition_number=singular[0]/singular[-1]))
    return rows


def load_data(cells=120):
    if not isinstance(cells, (int, np.integer)) or cells < 2 or cells % 2:
        raise ValueError("an even cell count >= 2 is required for the reflection basis")
    rows = []
    for geometry in ["lumped", "local"]:
        d = Device(D=.1, gamma=.2, kappa_left=1., kappa_right=1., a=2.,
                   cells=cells, geometry=geometry)
        baseline = d.decay_rates(count=cells)
        Y, _ = d.response(0.)
        c, K, B, _ = d.stiffness()
        # Reflection-even basis; only this sector couples to uniform charge drive.
        even = np.zeros((cells, cells//2))
        for i in range(cells//2):
            even[i, i] = even[cells-1-i, i] = 1/np.sqrt(2)
        Ke, be = even.T@K@even, even.T@B[:, 0]
        for R in np.logspace(-3, 3, 70):
            rates = d.decay_rates(R, count=cells)
            loaded_even = Ke+load_feedback(R, d.G)*np.outer(be, be)
            even_rate = eigh(loaded_even/c, subset_by_index=[0, 0], eigvals_only=True)[0]
            product = np.exp(np.sum(np.log(rates/baseline)))
            electrical = 1+load_feedback(R, d.G)*(Y[0, 0].real-d.G)
            rows.append(dict(geometry=geometry, RG=R*d.G, lowest_rate=rates[0],
                             lowest_even_rate=even_rate,
                             node_rate=node_rate(d, R), pole_product=product,
                             electrical_prediction=electrical))
    return rows


def port_data():
    rows = []
    for ratio in [1., .3, .1]:
        d = Device(kappa_right=.3*ratio)
        for w in np.logspace(-3, 3, 160):
            Y, _ = d.response(w)
            common, differential = Y[0, 1]+Y[0, 2], Y[0, 1]-Y[0, 2]
            rows.append(dict(contact_ratio=ratio, omega=w, common_abs=abs(common),
                             differential_abs=abs(differential),
                             differential_over_common=abs(differential/common)))
    return rows


def elimination_data():
    G = np.array([[100001., -10000.], [-10000., 1000.]])
    schur = G[0, 0]-G[0, 1]*G[1, 0]/G[1, 1]
    capacity = 1+G[0, 1]*G[1, 0]/G[1, 1]**2
    return dict(generator_per_second=G.tolist(), exact_rates=np.linalg.eigvalsh(G).tolist(),
                schur_rate=schur, capacity=capacity, corrected_rate=schur/capacity)


def main():
    for folder in ["figures", "data", "tables"]:
        (ROOT/folder).mkdir(exist_ok=True)
    samples, fits = contact_fits()
    information, loads, ports = length_information(), load_data(), port_data()
    for name, rows in [("contact_misspecification.csv", samples),
                       ("contact_misspecification_fits.csv", fits),
                       ("length_identifiability.csv", information),
                       ("geometry_load_comparison.csv", loads),
                       ("asymmetric_port_response.csv", ports)]:
        write_csv(name, rows)
    (ROOT/"data/dynamic_elimination.json").write_text(json.dumps(elimination_data(), indent=2)+"\n")
    table = [r"\begin{tabular}{lrrrr}", r"\toprule",
             r"Length range & $D_s$ & $\Gamma$ & $\kappa$ & RMS (\%) \\", r"\midrule"]
    for r in fits:
        table.append(f"${r['length_min']:g}$--${r['length_max']:g}$ & {r['D']:.3f} & "
                     f"{r['Gamma']:.3f} & {r['kappa']:.3f} & {100*r['rms_log_mismatch']:.3f} " + r"\\")
    table += [r"\bottomrule", r"\end{tabular}"]
    (ROOT/"tables/contact_misspecification.tex").write_text("\n".join(table)+"\n")
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    colors = ["#007c78", "#b83552", "#375b9a"]
    fig, ax = plt.subplots(1, 2, figsize=(7., 2.7), constrained_layout=True)
    broad = samples[:18]
    ax[0].loglog([r['length'] for r in broad], [r['asymmetric_rate'] for r in broad],
                 'o', markersize=3, color=colors[0], label="unequal contacts, exact")
    dense = np.geomspace(.05, 30., 250)
    p = fits[0]
    ax[0].loglog(dense, [slowest_pole(L, p['D'], p['Gamma'], p['kappa']) for L in dense],
                 '--', color=colors[1], label="equal-contact fit")
    ax[0].set(xlabel="Length (reference units)", ylabel="Lowest decay rate",
              title="(a) A good fit with biased diffusion")
    ax[0].legend(fontsize=7)
    ax[1].axhspan(-2, 2, color="#eeeeee", label="2% comparison scale")
    for i, r in enumerate(fits):
        selected = samples[18*i:18*(i+1)]
        ax[1].semilogx([v['length'] for v in selected], [100*v['log_residual'] for v in selected],
                       color=colors[i], label=f"L = {r['length_min']:g}-{r['length_max']:g}")
    ax[1].set(xlabel="Length (reference units)", ylabel="Log-rate residual (%)",
              title="(b) Contact-model mismatch")
    ax[1].legend(fontsize=7)
    fig.savefig(ROOT/"figures/fig_device_inference.pdf", metadata=METADATA)
    plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(7., 2.6), constrained_layout=True)
    for i, ratio in enumerate([1., .3, .1]):
        selected = [r for r in ports if r['contact_ratio'] == ratio]
        w = [r['omega'] for r in selected]
        label = rf"$\kappa_R/\kappa_L={ratio:g}$"
        ax[0].loglog(w, [r['common_abs'] for r in selected], color=colors[i], label=label)
        # Exact symmetry is shown as zero, not as solver roundoff.
        values = np.zeros(len(w)) if ratio == 1 else [r['differential_over_common'] for r in selected]
        ax[1].semilogx(w, values, color=colors[i], label=label)
    ax[0].set(xlabel=r"Frequency $\omega/\Gamma$", ylabel="Common inverse response",
              title="(a) Two-contact injection")
    ax[1].set(xlabel=r"Frequency $\omega/\Gamma$", ylabel="Differential / common magnitude",
              title="(b) Asymmetry lifts cancellation")
    for axis in ax:
        axis.legend(fontsize=7, loc="lower right")
    ax[0].legend(fontsize=7, loc="upper right")
    fig.savefig(ROOT/"figures/fig_reciprocal_ports.pdf", metadata=METADATA)
    plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(7., 2.7), constrained_layout=True)
    for i, geometry in enumerate(["lumped", "local"]):
        selected = [r for r in loads if r['geometry'] == geometry]
        x = [r['RG'] for r in selected]
        ax[0].semilogx(x, [r['lowest_rate'] for r in selected], color=colors[i], label=geometry)
        ax[1].semilogx(x, [r['electrical_prediction'] for r in selected], color=colors[i], label=geometry+", electrical")
        ax[1].semilogx(x[::7], [r['pole_product'] for r in selected][::7], 'o',
                       markersize=3, markerfacecolor="none", color=colors[i], label=geometry+", poles")
        if i == 0:
            ax[0].semilogx(x, [r['lowest_even_rate'] for r in selected], ':', color=colors[0],
                           label="lumped even mode")
            ax[0].semilogx(x, [r['node_rate'] for r in selected], '--', color=colors[2], label="uniform node")
    ax[0].set(xlabel=r"Electrical load $RG$", ylabel="Lowest rate (reference units)",
              title="(a) Geometry-dependent decay")
    ax[1].set(xlabel=r"Electrical load $RG$", ylabel="Loaded / short-circuit pole product",
              title="(b) Electrical consistency identity")
    for axis in ax:
        axis.legend(fontsize=7)
    fig.savefig(ROOT/"figures/fig_device_loading.pdf", metadata=METADATA)
    plt.close(fig)
    print("Generated robustness data, contact table, and three device figures.")


if __name__ == "__main__":
    main()
