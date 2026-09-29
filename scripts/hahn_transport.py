#!/usr/bin/env python3
"""Hahn-inspired helical Rashba spin alignment with a resolved thermal bath.

All rates are angular frequencies [1/s]. omega_c sets the bath cutoff;
theta = kB*T/(hbar*omega_c). Examples set omega_c=1 as a reference unit.
"""
from dataclasses import dataclass
import json
import numpy as np
from scipy.linalg import eigh
from reciprocal_device import ROOT, METADATA, plt, write_csv

SX = np.array([[0., 1.], [1., 0.]], complex)
SY = np.array([[0., -1j], [1j, 0.]], complex)
SZ = np.diag([1., -1.]).astype(complex)
IDENTITY = np.eye(2, dtype=complex)
SP = (SX+1j*SY)/2
SM = SP.conj().T


def flight_response(gamma, transit, frequency):
    """Exit/target response for fixed-rate modulation with exp(-i*frequency*t).

    gamma >= 0 and frequency are in 1/s; transit >= 0 is in seconds.
    The incident spin is unpolarized and the flight time is deterministic.
    """
    frequency = np.asarray(frequency, dtype=float)
    if not np.isfinite([gamma, transit]).all() or min(gamma, transit) < 0 or not np.isfinite(frequency).all():
        raise ValueError("finite nonnegative rate/time and finite frequencies required")
    if gamma == 0:
        return np.zeros_like(frequency, dtype=complex)
    q = gamma-1j*frequency
    return -gamma*np.expm1(-q*transit)/q


@dataclass(frozen=True)
class Helix:
    """Selected trajectory: Omega=chi*flow*omega, a=flow*zeta*omega.

    omega>0 [1/s], zeta is a signed dimensionless mean SOC, pitch_ratio=radius/b>0,
    g>=0 [1/s] fluctuating SOC amplitude. chi/flow are +/-1.
    Setting all SOC to zero requires both zeta=0 and g=0.
    """
    omega: float = 1.
    zeta: float = .15
    pitch_ratio: float = .7
    g: float = .08
    omega_c: float = 1.
    theta: float = .2
    chi: int = 1
    flow: int = 1

    def __post_init__(self):
        values = [self.omega, self.zeta, self.pitch_ratio, self.g, self.omega_c, self.theta]
        if not np.isfinite(values).all() or min(self.omega, self.pitch_ratio, self.omega_c, self.theta) <= 0:
            raise ValueError("positive finite frequency, pitch ratio, cutoff and temperature required")
        if self.g < 0 or self.chi not in (-1, 1) or self.flow not in (-1, 1):
            raise ValueError("nonnegative fluctuation amplitude and signed chirality/flow required")

    @property
    def rotation(self):
        return self.chi*self.flow*self.omega

    @property
    def vector(self):
        a = self.flow*self.zeta*self.omega
        return np.array([0., -a, self.chi*self.pitch_ratio*a-self.rotation])

    @property
    def axis(self):
        return self.vector/np.linalg.norm(self.vector)

    @property
    def coupling(self):
        return (-SY+self.chi*self.pitch_ratio*SZ)/np.sqrt(1+self.pitch_ratio**2)

    @property
    def hamiltonian(self):
        """H_rot/hbar, a 2-by-2 Hermitian rate matrix."""
        return (self.vector[1]*SY+self.vector[2]*SZ)/2

    def unitary(self, time):
        angle = self.rotation*time
        return np.diag([np.exp(-.5j*angle), np.exp(.5j*angle)])

    def spectrum(self, frequency):
        """Unsymmetrized bath spectrum S_X [s], with thermal KMS balance."""
        if not np.isfinite(frequency):
            raise ValueError("finite spectral frequency required")
        x = abs(frequency)/self.omega_c
        if x == 0:
            return self.theta/self.omega_c
        positive = x*np.exp(-x)/(-np.expm1(-x/self.theta))/self.omega_c
        return positive if frequency > 0 else positive*np.exp(-x/self.theta)

    def rates(self, bath="covariant"):
        """Return downhill/uphill rates and energy eigenvectors (lower,upper).

        The covariant Rashba-amplitude bath is static in the rotating frame.
        A lab-fixed sigma_x bath instead has m=+/-1 Floquet sidebands.
        """
        values, vectors = eigh(self.hamiltonian)
        gap = values[1]-values[0]
        if bath == "covariant":
            components = [(0, self.coupling)]
        elif bath == "lab_fixed":
            components = [(1, SP), (-1, SM)]
        else:
            raise ValueError("bath must be covariant or lab_fixed")
        down, up = 0., 0.
        for order, operator in components:
            matrix = vectors.conj().T@operator@vectors
            down += abs(matrix[0, 1])**2*self.spectrum(gap-order*self.rotation)
            up += abs(matrix[1, 0])**2*self.spectrum(-gap-order*self.rotation)
        return self.g**2*down/4, self.g**2*up/4, vectors

    def polarization(self, time, bath="covariant"):
        time = np.asarray(time, dtype=float)
        if not np.isfinite(time).all() or (time < 0).any():
            raise ValueError("finite nonnegative evolution times required")
        down, up, _ = self.rates(bath)
        total = down+up
        if total == 0:
            return np.zeros_like(time)
        target = self.axis[2]*(up-down)/total
        return target*(-np.expm1(-total*time))

    def target(self, bath="covariant"):
        down, up, _ = self.rates(bath)
        return np.nan if down+up == 0 else self.axis[2]*(up-down)/(down+up)

    def jumps(self):
        """Rotating-frame secular GKSL jumps [s^-1/2], including dephasing."""
        down, up, vectors = self.rates()
        lower, upper = vectors[:, 0], vectors[:, 1]
        diagonal = np.diag(vectors.conj().T@self.coupling@vectors).real
        dephase = self.g*np.sqrt(self.spectrum(0))/2 * (vectors@np.diag(diagonal)@vectors.conj().T)
        return [np.sqrt(down)*np.outer(lower, upper.conj()),
                np.sqrt(up)*np.outer(upper, lower.conj()), dephase]


def dissipator(state, jumps):
    result = np.zeros((2, 2), complex)
    for jump in jumps:
        square = jump.conj().T@jump
        result += jump@state@jump.conj().T-(square@state+state@square)/2
    return result


def example_data():
    base = Helix()
    down, up, _ = base.rates()
    buildup = [dict(time_over_T1=x, Pz=float(base.polarization(x/(down+up))),
                    target=base.target()) for x in np.linspace(0, 6, 241)]
    soc = []
    for scale in np.linspace(0, 1, 201):
        model = Helix(zeta=.15*scale, g=.08*scale)
        d, u, _ = model.rates()
        for transit in [100., 1000., 10000.]:
            soc.append(dict(SOC_scale=scale, transit_omega_c=transit,
                            Gamma1_over_omega_c=d+u, Pz=float(model.polarization(transit))))
    baths = []
    for omega in np.geomspace(.02, 4., 180):
        model = Helix(omega=omega, g=.04*omega)
        for bath in ["covariant", "lab_fixed"]:
            d, u, _ = model.rates(bath)
            baths.append(dict(omega_over_cutoff=omega, bath=bath, target=model.target(bath),
                              Gamma1_over_omega_c=d+u,
                              Pz_ten_turns=float(model.polarization(20*np.pi/omega, bath))))
    summary = dict(reference_units="omega_c=1; theta=kB*T/(hbar*omega_c)=0.2",
                   base=base.__dict__, base_gap=float(np.linalg.norm(base.vector)),
                   base_Gamma1=down+up, base_target=base.target(),
                   SOC_scan=dict(zeta_max=.15, g_max=.08, transit_times=[100., 1000., 10000.]),
                   bath_scan=dict(g_over_omega=.04, helix_turns=10),
                   status="synthetic mechanism calculations, not a fit to molecular data")
    return buildup, soc, baths, summary


def main():
    for folder in ("data", "figures"):
        (ROOT/folder).mkdir(exist_ok=True)
    buildup, soc, baths, summary = example_data()
    write_csv("hahn_alignment.csv", buildup)
    write_csv("hahn_soc_transit.csv", soc)
    write_csv("hahn_bath_comparison.csv", baths)
    (ROOT/"data/hahn_parameters.json").write_text(json.dumps(summary, indent=2)+"\n")
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(7.1, 3.0), layout="constrained")
    ax = fig.add_subplot(121, projection="3d")
    phi = np.linspace(0, 4*np.pi, 350)
    ax.plot(.7*np.cos(phi), .7*np.sin(phi), phi/(2*np.pi), color="#287c83", lw=2)
    ax.quiver(0, 0, .1, 0, 0, 1.6, color="#ad3d38", linewidth=2, arrow_length_ratio=.12)
    ax.text(0, 0, 1.9, r"$\Omega$, spin axis", fontsize=9)
    ax.set(xlabel="x / b", ylabel="y / b", zlabel="z / pitch", title="(a) Directed helical motion")
    ax.set_xticks([-.7, .7]); ax.set_yticks([-.7, .7]); ax.set_zticks([0, 1, 2])
    ax.view_init(elev=16, azim=-55)
    ax = fig.add_subplot(122)
    ax.plot([r["time_over_T1"] for r in buildup], [r["Pz"] for r in buildup], color="#287c83", lw=2)
    ax.axhline(summary["base_target"], color="#ad3d38", ls="--", label="Alignment target")
    ax.set(xlabel=r"Interaction time $t/T_{1,\mathrm{rot}}$", ylabel=r"Physical-spin polarization $P_z$",
           title="(b) Relaxation builds alignment", ylim=(0, 1.05))
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(ROOT/"figures/fig_hahn_alignment.pdf", metadata=METADATA)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.8), layout="constrained")
    for transit, color in [(100., "#287c83"), (1000., "#ad3d38"), (10000., "#58559a")]:
        rows = [r for r in soc if r["transit_omega_c"] == transit]
        axes[0].plot([r["SOC_scale"] for r in rows], [r["Pz"] for r in rows], color=color,
                     label=rf"$\omega_c\tau={transit:g}$")
    rows = [r for r in soc if r["transit_omega_c"] == 100.]
    axes[1].plot([r["SOC_scale"] for r in rows], [r["Gamma1_over_omega_c"] for r in rows], color="#287c83")
    axes[0].set(xlabel=r"SOC scale $\lambda$", ylabel=r"Transmitted polarization $P_z^{\rm out}$",
                title="(a) Finite-transit spin selectivity", ylim=(0, 1.05))
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set(xlabel=r"SOC scale $\lambda$", ylabel=r"Alignment rate $1/(\omega_c T_{1,\mathrm{rot}})$",
                title="(b) SOC-controlled relaxation")
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.savefig(ROOT/"figures/fig_hahn_transit.pdf", metadata=METADATA)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.8), layout="constrained")
    for bath, label, color in [("covariant", "Rashba-amplitude bath", "#287c83"),
                                ("lab_fixed", "Lab-fixed transverse bath", "#ad3d38")]:
        rows = [r for r in baths if r["bath"] == bath]
        omega = [r["omega_over_cutoff"] for r in rows]
        axes[0].semilogx(omega, [r["target"] for r in rows], label=label, color=color)
        axes[1].semilogx(omega, [r["Pz_ten_turns"] for r in rows], color=color)
    axes[0].set(xlabel=r"Traversal rotation $|\Omega|/\omega_c$", ylabel=r"Stationary $P_z$",
                title="(a) Bath coupling sets the target")
    axes[0].legend(frameon=False, fontsize=7.5)
    axes[1].set(xlabel=r"Traversal rotation $|\Omega|/\omega_c$", ylabel=r"$P_z^{\rm out}$ after ten turns",
                title="(b) Finite-time alignment")
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.savefig(ROOT/"figures/fig_hahn_baths.pdf", metadata=METADATA)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
