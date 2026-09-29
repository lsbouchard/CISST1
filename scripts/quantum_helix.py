#!/usr/bin/env python3
"""Static helical Rashba completion and quantitative trajectory validity tests.

Quantum matrices are H/hbar [1/s], p_phi/hbar and Q/hbar (dimensionless).
The Fourier truncation is an operator check, not a two-lead scattering model.
"""
from dataclasses import replace
import json
import numpy as np
from hahn_transport import Helix, SY, SZ
from generate_figures import HBAR, KB, GAMMA_E_ABS, N_TURN
from reciprocal_device import ROOT, METADATA, plt, write_csv

ELECTRON_MASS = 9.1093837139e-31  # CODATA 2022 [kg]; illustrative m_* only.


def fourier_operators(cutoff, hbar_over_inertia, zeta, pitch_ratio, chi=1):
    """Galerkin matrices on |l,s>, l=-cutoff,...,cutoff, s=up/down.

    hbar_over_inertia>0 [1/s], signed zeta, pitch_ratio>0, chi=+/-1.
    Returns rate Hamiltonian, dimensionless screw generator and physical Sz/hbar.
    """
    if not isinstance(cutoff, (int, np.integer)) or cutoff < 1:
        raise ValueError("positive integer Fourier cutoff required")
    if not np.isfinite([hbar_over_inertia, zeta, pitch_ratio]).all() or min(hbar_over_inertia, pitch_ratio) <= 0 or chi not in (-1, 1):
        raise ValueError("finite SOC, positive inertia frequency/geometry, signed chirality required")
    ell = np.arange(-cutoff, cutoff+1)
    p = np.kron(np.diag(ell), np.eye(2))
    spin = np.kron(np.eye(len(ell)), SZ/2)
    A = np.kron(np.eye(len(ell)), chi*pitch_ratio*SZ)
    # -sigma_phi connects |l,down> to |l-1,up> with amplitude +i.
    for i in range(1, len(ell)):
        A[2*(i-1), 2*i+1] = 1j
        A[2*i+1, 2*(i-1)] = -1j
    H = hbar_over_inertia*(p@p/2+chi*zeta*(A@p+p@A)/4)
    return H, p+spin, spin


def sector_hamiltonian(k, hbar_over_inertia, zeta, pitch_ratio, chi=1):
    """Exact 2x2 fixed-Q/hbar=k bulk Hamiltonian H_k/hbar [1/s]."""
    if not np.isfinite([k, hbar_over_inertia, zeta, pitch_ratio]).all() or min(hbar_over_inertia, pitch_ratio) <= 0 or chi not in (-1, 1):
        raise ValueError("finite k/SOC, positive inertia frequency/geometry, signed chirality required")
    scalar = hbar_over_inertia*((k*k+.25)/2-zeta*pitch_ratio/4)
    omega = hbar_over_inertia*k
    return scalar*np.eye(2)+( -chi*zeta*omega*SY+(pitch_ratio*zeta-1)*omega*SZ)/2


def sector_polarization(model, time, initial_spin):
    """Physical Pz for initial +/-z in one nonzero-Q sector of the GKSL model."""
    if initial_spin not in (-1, 1):
        raise ValueError("initial spin must be +1 or -1")
    time = np.asarray(time, dtype=float)
    if not np.isfinite(time).all() or (time < 0).any():
        raise ValueError("finite nonnegative times required")
    down, up, _ = model.rates()
    gap = np.linalg.norm(model.vector)
    nz = model.axis[2]
    m = np.array([0., -1., model.chi*model.pitch_ratio])/np.sqrt(1+model.pitch_ratio**2)
    gamma_phi = model.g**2*np.dot(m, model.axis)**2*model.spectrum(0)/2
    gamma1 = down+up
    gamma2 = gamma1/2+gamma_phi
    thermal = np.tanh(gap/(2*model.theta*model.omega_c))
    return (-nz*thermal+nz*(initial_spin*nz+thermal)*np.exp(-gamma1*time)
            +initial_spin*(1-nz*nz)*np.exp(-gamma2*time)*np.cos(gap*time))


def common_momentum_polarization(ell, time, model=Helix()):
    """Equal mixture with |p_phi|/hbar=ell>0, ell!=1/2, at entry.

    The sign of p_phi follows model.rotation; Q/hbar=p_phi/hbar+/-1/2.
    Mean traversal rotation magnitude is held fixed as ell changes.
    The fluctuating Rashba coefficient, not g, is the same in both sectors.
    This returns spin at fixed time, not a spin-resolved transmission flux.
    """
    if not np.isfinite(ell) or ell <= 0 or ell == .5:
        raise ValueError("positive orbital action, excluding the zero-gap sector at ell=1/2, required")
    answer = 0.
    momentum = np.sign(model.rotation)*ell
    for spin in (-1, 1):
        factor = (momentum+spin/2)/momentum
        sector = replace(model, omega=model.omega*abs(factor), g=model.g*abs(factor),
                         flow=int(model.flow*np.sign(factor)))
        answer = answer+sector_polarization(sector, time, spin)/2
    return answer


def alignment_turns(fraction, rate_over_rotation):
    """Turns needed to attain 0<fraction<1 of |P_*| in constant flight."""
    rate_over_rotation = np.asarray(rate_over_rotation, dtype=float)
    if not np.isfinite(fraction) or not 0 < fraction < 1 or not np.isfinite(rate_over_rotation).all() or (rate_over_rotation <= 0).any():
        raise ValueError("fraction in (0,1) and finite positive rate ratios required")
    return -np.log1p(-fraction)/(2*np.pi*rate_over_rotation)


def example_data():
    base = Helix()
    down, up, _ = base.rates()
    ratio = (down+up)/base.omega
    turns = [dict(Gamma1_over_abs_Omega=q, alignment_fraction=f,
                  required_turns=float(alignment_turns(f, q)))
             for f in [.1, .5, .9] for q in np.geomspace(1e-5, .1, 181)]
    radius, pitch = 1e-9, 3.4e-9
    b = pitch/(2*np.pi)
    inertia = ELECTRON_MASS*(radius**2+b**2)
    recoil = HBAR**2/(2*inertia)
    scales = []
    for rate in np.geomspace(1e11, 1e16, 201):
        omega = 2*np.pi*rate/N_TURN
        scales.append(dict(assumed_step_rate_per_s=rate, abs_Omega_per_s=omega,
                           B_equivalent_T=omega/GAMMA_E_ABS,
                           orbital_action_over_hbar=inertia*omega/HBAR,
                           recoil_energy_J=recoil, spin_frame_gap_J=HBAR*omega,
                           kinematic_target_at_300K=np.tanh(HBAR*omega/(2*KB*300))))
    quantum = []
    for ell in [1., 2., 5., 10., 100., 1000.]:
        for time in [0., 20*np.pi, 1/(down+up), 4/(down+up)]:
            exact = float(common_momentum_polarization(ell, time))
            trajectory = float(base.polarization(time))
            quantum.append(dict(orbital_action_over_hbar=ell, time_omega_c=time,
                                common_momentum_Pz=exact, trajectory_Pz=trajectory,
                                absolute_difference=abs(exact-trajectory)))
    summary = dict(status="conditional bulk quantum completion; no reservoir transmission fit",
                   constants=dict(hbar_J_s=HBAR, kB_J_per_K=KB, gamma_e_abs_per_s_T=GAMMA_E_ABS,
                                  steps_per_turn=N_TURN, illustrative_mass_kg=ELECTRON_MASS),
                   illustrative_geometry=dict(radius_m=radius, pitch_m=pitch, inertia_kg_m2=inertia),
                   recoil_energy_J=recoil,
                   orbital_action_at_step_rate_1e13=inertia*(2*np.pi*1e13/N_TURN)/HBAR,
                   reference_Gamma1_over_abs_Omega=ratio,
                   reference_turns_per_T1=1/(2*np.pi*ratio),
                   reference_turns_to_90_percent=float(alignment_turns(.9, ratio)),
                   reference_ten_turn_Pz=float(base.polarization(20*np.pi/base.omega)))
    return turns, scales, quantum, summary


def main():
    turns, scales, quantum, summary = example_data()
    write_csv("alignment_turns.csv", turns)
    write_csv("orbital_recoil_scales.csv", scales)
    write_csv("quantum_helix_matching.csv", quantum)
    (ROOT/"data/quantum_helix_parameters.json").write_text(json.dumps(summary, indent=2)+"\n")
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.9), layout="constrained")
    for fraction, color in [(.1, "#287c83"), (.5, "#ad3d38"), (.9, "#58559a")]:
        rows = [r for r in turns if r["alignment_fraction"] == fraction]
        axes[0].loglog([r["Gamma1_over_abs_Omega"] for r in rows],
                       [r["required_turns"] for r in rows], color=color, label=f"{100*fraction:g}% of target")
    axes[0].axvline(summary["reference_Gamma1_over_abs_Omega"], color=".5", ls=":")
    axes[0].axhline(10, color=".5", ls="--")
    axes[0].set(xlabel=r"Alignment rate / traversal rate $\Gamma_1/|\Omega|$",
                ylabel="Required active helical turns", title="(a) Time required for alignment")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].loglog([r["assumed_step_rate_per_s"] for r in scales],
                   [r["orbital_action_over_hbar"] for r in scales], color="#287c83")
    axes[1].axhline(1, color=".5", ls="--")
    axes[1].axvline(1e13, color=".5", ls=":")
    axes[1].set(xlabel=r"Assumed step rate $R_{\rm hop}$ (s$^{-1}$)",
                ylabel=r"Orbital action $I_h|\Omega|/\hbar$",
                title=r"(b) Semiclassical consistency ($m_*=m_e$)")
    fig.savefig(ROOT/"figures/fig_alignment_validity.pdf", metadata=METADATA)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
