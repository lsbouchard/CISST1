#!/usr/bin/env python3
"""Conserving finite helical junction with a thermal Rashba-amplitude bath.

Reference units: energy t, time hbar/t, particle current t/hbar. All matrices
act on C^(2*sites), ordered (site, physical spin up/down). No Barnett Zeeman
term or spin-relaxation target is added to the static laboratory Hamiltonian.
"""
from dataclasses import dataclass
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import expit
from hahn_transport import SX, SY, SZ, IDENTITY


def adjoint(a):
    return np.swapaxes(a.conj(), -1, -2)


def hermitian(a):
    return (a+adjoint(a))/2


@dataclass(frozen=True)
class Junction:
    """Static lattice realization of the symmetrized helical Hamiltonian.

    sites>=3, angular step>0, t>0, signed mean zeta and chirality +/-1.
    The coordinate is axial angle u=z/b; actual azimuth is phi=chi*u.
    Fluctuations zeta -> zeta+eta*X use the exact derivative of this H.
    Lead hopping and contact hoppings are positive energies in the same units.
    """
    sites: int = 9
    step: float = np.pi/4
    t: float = 1.
    zeta: float = .15
    pitch_ratio: float = .7
    chi: int = 1
    eta: float = .08
    lead_hopping: float = 2.
    contact_left: float = 1.
    contact_right: float = 1.

    def __post_init__(self):
        values = [self.step, self.t, self.zeta, self.pitch_ratio, self.eta,
                  self.lead_hopping, self.contact_left, self.contact_right]
        if not isinstance(self.sites, (int, np.integer)) or self.sites < 3:
            raise ValueError("at least three integer sites required")
        if not np.isfinite(values).all() or min(self.step, self.t, self.pitch_ratio,
                self.lead_hopping, self.contact_left, self.contact_right) <= 0 or self.eta < 0 or self.chi not in (-1, 1):
            raise ValueError("finite SOC, positive geometry/hoppings, nonnegative noise, signed chirality required")

    @property
    def turns(self):
        return (self.sites-1)*self.step/(2*np.pi)

    def operators(self):
        dimension = 2*self.sites
        H = np.eye(dimension, dtype=complex)*(2*self.t)
        derivative = np.zeros_like(H)
        phi = self.chi*(np.arange(self.sites)-(self.sites-1)/2)*self.step
        A = np.sin(phi)[:, None, None]*SX-np.cos(phi)[:, None, None]*SY+self.chi*self.pitch_ratio*SZ
        for j in range(self.sites-1):
            a, b = slice(2*j, 2*j+2), slice(2*j+2, 2*j+4)
            vertex = -.5j*self.t*self.step*(A[j]+A[j+1])/2
            derivative[a, b], derivative[b, a] = vertex, vertex.conj().T
            H[a, b] = -self.t*IDENTITY+self.zeta*vertex
            H[b, a] = H[a, b].conj().T
        return H, self.eta*derivative

    def embeddings(self, energies):
        """Retarded lead self-energies and PSD broadenings, shape (2,NE,D,D)."""
        energies = np.asarray(energies, dtype=float)
        if not np.isfinite(energies).all():
            raise ValueError("finite energies required")
        x, t = energies-2*self.t, self.lead_hopping
        root = np.sqrt(np.maximum(x*x-4*t*t, 0.))
        g = (x-np.sign(x)*root)/(2*t*t)-1j*np.sqrt(np.maximum(4*t*t-x*x, 0.))/(2*t*t)
        out = np.zeros((2, len(energies), 2*self.sites, 2*self.sites), complex)
        for port, hopping, site in [(0, self.contact_left, 0), (1, self.contact_right, self.sites-1)]:
            block = slice(2*site, 2*site+2)
            out[port, :, block, block] = hopping*hopping*g[:, None, None]*IDENTITY
        return out, hermitian(1j*(out-adjoint(out)))


@dataclass(frozen=True)
class ThermalBath:
    """X=sum c_j(b_j+b_j^dagger), energies hbar*omega_j and c_j^2 weights.

    Temperature means kB*T in energy units. Positive weights are dimensionless.
    The phonons are maintained at this temperature; their heat flow is measured.
    """
    energies: tuple = (.4,)
    weights: tuple = (1.,)
    temperature: float = .15

    def __post_init__(self):
        e, w = np.asarray(self.energies), np.asarray(self.weights)
        if e.ndim != 1 or len(e) == 0 or e.shape != w.shape or not np.isfinite(e).all() or not np.isfinite(w).all() or (e <= 0).any() or (w <= 0).any() or not np.isfinite(self.temperature) or self.temperature < 0:
            raise ValueError("positive finite phonon energies/weights and nonnegative temperature required")

    @property
    def occupations(self):
        if self.temperature == 0:
            return np.zeros(len(self.energies))
        x = np.asarray(self.energies)/self.temperature
        return np.exp(-x)/(-np.expm1(-x))

    @classmethod
    def ohmic(cls, spacing=.1, cutoff=.5, count=30, temperature=.15):
        """Quadrature for S_X(+w)=w*exp(-w/wc)/(wc^2*(1-exp(-w/T)))."""
        if not np.isfinite([spacing, cutoff]).all() or min(spacing, cutoff) <= 0 or not isinstance(count, (int, np.integer)) or count < 1:
            raise ValueError("positive quadrature spacing/cutoff and integer count required")
        e = spacing*np.arange(1, count+1)
        weights = spacing*e*np.exp(-e/cutoff)/(2*np.pi*cutoff**2)
        return cls(tuple(e), tuple(weights), temperature)


def energy_grid(low=-5., high=9., spacing=.025):
    if not np.isfinite([low, high, spacing]).all() or spacing <= 0 or high <= low:
        raise ValueError("increasing finite energy interval and positive spacing required")
    count = int(round((high-low)/spacing))
    if not np.isclose(count*spacing, high-low, rtol=1e-12, atol=1e-12):
        raise ValueError("energy interval must be an integer number of spacings")
    return low+spacing*np.arange(count+1)


def shifted(a, offset):
    """a(E+offset*dE), with zero spectral weight beyond the numerical window."""
    out = np.zeros_like(a)
    if offset == 0:
        return a.copy()
    if abs(offset) >= len(a):
        return out
    if offset > 0:
        out[:-offset] = a[offset:]
    else:
        out[-offset:] = a[:offset]
    return out


def hermitian_convolution(a, kernel):
    """Full-matrix real-kernel convolution, packing only Hermitian redundancy.

    Every spatial/spin off-diagonal entry is retained. Contiguous energy traces
    avoid slow strided FFTs; the lower triangle is fixed by Hermiticity.
    """
    dimension = a.shape[-1]
    row, column = np.triu_indices(dimension)
    packed = np.ascontiguousarray(a[:, row, column].T)
    convolved = fftconvolve(packed, kernel[None, :], mode="same", axes=-1).T
    out = np.empty_like(a)
    out[:, row, column] = convolved
    out[:, column, row] = convolved.conj()
    diagonal = np.arange(dimension)
    out[:, diagonal, diagonal] = out[:, diagonal, diagonal].real
    return out


def principal_value(broadening):
    """Exact PV of piecewise-linear Gamma, zero one grid step outside its window."""
    size = len(broadening)
    index = np.arange(-size+1, size, dtype=float)
    def xlogx(x):
        out = np.zeros_like(x)
        nonzero = x != 0
        out[nonzero] = x[nonzero]*np.log(abs(x[nonzero]))
        return out
    kernel = (xlogx(index+1)-2*xlogx(index)+xlogx(index-1))/(2*np.pi)
    return hermitian_convolution(broadening, kernel)


def phonon_self_energies(occupied, empty, vertex, bath, spacing):
    """Fock SCBA: positive in/out self-energies and causal retarded matrix."""
    offsets = np.rint(np.asarray(bath.energies)/spacing).astype(int)
    if (offsets < 1).any() or not np.allclose(offsets*spacing, bath.energies, atol=1e-12, rtol=1e-12):
        raise ValueError("every phonon energy must be an integer multiple of dE")
    center = int(max(offsets))
    kernel = np.zeros(2*center+1)
    for energy, weight, n in zip(bath.energies, bath.weights, bath.occupations):
        offset = int(round(energy/spacing))
        kernel[center-offset] += weight*(n+1)
        kernel[center+offset] += weight*n
    inside = hermitian_convolution(occupied, kernel)
    outside = hermitian_convolution(empty, kernel[::-1])
    inside = hermitian(vertex@inside@vertex)
    outside = hermitian(vertex@outside@vertex)
    broadening = inside+outside
    return inside, outside, principal_value(broadening)-.5j*broadening


def hartree_shift(density, vertex, bath):
    """Static mean displacement, -2*sum(c_j^2/w_j)*Tr(M*rho)*M [energy]."""
    coefficient = np.sum(np.asarray(bath.weights)/np.asarray(bath.energies))
    return -2*coefficient*np.trace(vertex@density).real*vertex


def phonon_heat_events(occupied, empty, vertex, bath, spacing, energies):
    """Heat to phonons from separately counted emission minus absorption events."""
    heat = 0.
    for energy, weight, n in zip(bath.energies, bath.weights, bath.occupations):
        offset = int(round(energy/spacing))
        emission = np.einsum('eij,eji->e', vertex@shifted(empty, -offset)@vertex, occupied).real
        absorption = np.einsum('eij,eji->e', vertex@shifted(empty, offset)@vertex, occupied).real
        heat += energy*weight*np.trapezoid((n+1)*emission-n*absorption, energies)/(2*np.pi)
    return float(heat)


@dataclass
class Solution:
    model: Junction
    bath: ThermalBath
    energies: np.ndarray
    chemical_potentials: tuple
    lead_temperature: float
    reservoir_kind: str
    retarded: np.ndarray
    occupied: np.ndarray
    empty: np.ndarray
    lead_retarded: np.ndarray
    lead_in: np.ndarray
    lead_out: np.ndarray
    phonon_retarded: np.ndarray
    phonon_in: np.ndarray
    phonon_out: np.ndarray
    hartree: np.ndarray
    iterations: int
    residual: float

    @property
    def density(self):
        return hermitian(np.trapezoid(self.occupied, self.energies, axis=0)/(2*np.pi))

    def observables(self):
        """Lead currents are INTO the channel; spin-rate uses sigma_z, not Sz.

        Net drain ratio and gross extraction polarization are distinct outputs.
        Convert sigma_z rate to spin angular momentum current by hbar/2.
        """
        spin = np.kron(np.eye(self.model.sites), SZ)
        incoming = np.einsum('peij,eji->pe', self.lead_in, self.empty).real
        outgoing = np.einsum('peij,eji->pe', self.lead_out, self.occupied).real
        spin_in = np.einsum('i,peij,eji->pe', np.diag(spin), self.lead_in, self.empty).real
        spin_out = np.einsum('i,peij,eji->pe', np.diag(spin), self.lead_out, self.occupied).real
        integrate = lambda a: np.trapezoid(a, self.energies, axis=-1)/(2*np.pi)
        particle, spin_rate = integrate(incoming-outgoing), integrate(spin_in-spin_out)
        gross_in, gross_out = integrate(incoming), integrate(outgoing)
        gross_spin = integrate(spin_out)
        electronic_energy = integrate((incoming-outgoing)*self.energies)
        collision = np.einsum('eij,eji->e', self.phonon_in, self.empty).real-np.einsum('eij,eji->e', self.phonon_out, self.occupied).real
        bath_number = float(integrate(collision))
        bath_power_into_electrons = float(integrate(collision*self.energies))
        H, vertex = self.model.operators()
        rho = self.density
        GA = adjoint(self.retarded)
        bath_matrix = -1j*(self.phonon_retarded@self.occupied-self.occupied@adjoint(self.phonon_retarded)
                          +self.phonon_in@GA-self.retarded@self.phonon_in)
        bath_torque = np.trace((spin/2)@np.trapezoid(bath_matrix, self.energies, axis=0)).real/(2*np.pi)
        hamiltonian_torque = np.trace(1j*((H+self.hartree)@(spin/2)-(spin/2)@(H+self.hartree))@rho).real
        bath_variance = float(np.sum(np.asarray(self.bath.weights)*(2*self.bath.occupations+1)))
        ph_norm = max(np.linalg.norm(s, 2) for s in self.phonon_retarded)/self.model.t
        current_scale = max(np.max(abs(particle)), np.max(gross_in+gross_out), 1e-15)
        extraction = np.divide(gross_spin, gross_out, out=np.zeros(2), where=gross_out > 1e-15)
        ratio = np.divide(spin_rate, particle, out=np.full(2, np.nan), where=abs(particle) > 1e-12)
        accumulated = np.trace(spin@rho).real/np.trace(rho).real if np.trace(rho).real > 1e-15 else 0.
        spectral = hermitian(1j*(self.retarded-adjoint(self.retarded)))
        sumrule = hermitian(np.trapezoid(spectral, self.energies, axis=0)/(2*np.pi))
        heat_events = phonon_heat_events(self.occupied, self.empty, vertex, self.bath,
                                         self.energies[1]-self.energies[0], self.energies)
        density_eigenvalues = np.linalg.eigvalsh(rho)
        return dict(particle_into_left=float(particle[0]), particle_into_right=float(particle[1]),
                    spin_sigma_into_left=float(spin_rate[0]), spin_sigma_into_right=float(spin_rate[1]),
                    net_spin_per_particle_left=float(ratio[0]) if np.isfinite(ratio[0]) else None,
                    net_spin_per_particle_right=float(ratio[1]) if np.isfinite(ratio[1]) else None,
                    extraction_Pz_left=float(extraction[0]), extraction_Pz_right=float(extraction[1]),
                    gross_extraction_left=float(gross_out[0]), gross_extraction_right=float(gross_out[1]),
                    gross_injection_left=float(gross_in[0]), gross_injection_right=float(gross_in[1]),
                    accumulated_Pz=float(accumulated), electron_number=float(np.trace(rho).real),
                    charge_balance_relative=float(abs(particle.sum())/current_scale),
                    phonon_number_collision=float(bath_number),
                    energy_balance=float(electronic_energy.sum()+bath_power_into_electrons),
                    energy_into_left=float(electronic_energy[0]), energy_into_right=float(electronic_energy[1]),
                    heat_into_phonons=float(-bath_power_into_electrons),
                    heat_from_events=heat_events,
                    event_energy_balance=float(electronic_energy.sum()-heat_events),
                    spin_balance=float(spin_rate.sum()/2+hamiltonian_torque+bath_torque),
                    hamiltonian_spin_torque=float(hamiltonian_torque), bath_spin_torque=float(bath_torque),
                    phonon_self_energy_over_t=float(ph_norm), bath_variance=float(bath_variance),
                    noise_vertex_rms_over_t=float(np.linalg.norm(vertex, 2)*np.sqrt(bath_variance)/self.model.t),
                    spectral_sumrule_error=float(np.linalg.norm(sumrule-np.eye(len(H)), 2)),
                    occupied_empty_identity_error=float(np.max(abs(self.occupied+self.empty-spectral))),
                    density_min_eigenvalue=float(density_eigenvalues[0]),
                    density_max_eigenvalue=float(density_eigenvalues[-1]),
                    iterations=self.iterations, fixed_point_residual=self.residual)


def solve(model=Junction(), bath=None, energies=None, chemical_potentials=(1.4, .6),
          lead_temperature=None, tolerance=2e-9, max_iterations=500, mixing=.35, include_hartree=True,
          basis=None, lead_occupations=None):
    """Full spin-matrix Fock+Hartree SCBA with Fermi leads and thermal phonons.

    Convergence is checked against the unmixed self-energy fixed point. Energy
    window, grid, bath quadrature and approximation errors need separate tests.
    Optional spin-independent occupations (2,NE), between 0 and 1, permit a
    controlled dilute-beam comparison; those reservoirs are not Fermi baths.
    """
    bath = ThermalBath.ohmic() if bath is None else bath
    energies = energy_grid() if energies is None else np.asarray(energies, dtype=float)
    temperature = bath.temperature if lead_temperature is None else lead_temperature
    if energies.ndim != 1 or len(energies) < 5 or not np.isfinite(energies).all() or not np.isfinite([*chemical_potentials, temperature, tolerance, mixing]).all() or temperature < 0 or tolerance <= 0 or not 0 < mixing <= 1 or len(chemical_potentials) != 2 or not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
        raise ValueError("valid uniform grid, two finite potentials, and positive solver controls required")
    spacing = energies[1]-energies[0]
    if spacing <= 0 or not np.allclose(np.diff(energies), spacing, rtol=1e-11, atol=1e-12):
        raise ValueError("strictly increasing uniform energy grid required")
    H, M = model.operators()
    lead_R, gamma = model.embeddings(energies)
    if basis is not None:
        basis = np.asarray(basis, dtype=complex)
        if basis.shape != H.shape or not np.isfinite(basis).all() or not np.allclose(adjoint(basis)@basis, np.eye(len(H)), atol=1e-12, rtol=0):
            raise ValueError("basis must be a unitary 2*sites square matrix")
        H, M = adjoint(basis)@H@basis, adjoint(basis)@M@basis
        lead_R, gamma = adjoint(basis)@lead_R@basis, adjoint(basis)@gamma@basis
    f = np.array([expit((mu-energies)/temperature) if temperature > 0 else np.where(energies < mu, 1., np.where(energies > mu, 0., .5))
                  for mu in chemical_potentials])
    reservoir_kind = 'fermi'
    if lead_occupations is not None:
        f = np.asarray(lead_occupations, dtype=float)
        if f.shape != (2, len(energies)) or not np.isfinite(f).all() or (f < 0).any() or (f > 1).any():
            raise ValueError("spin-independent lead occupations must have shape (2,NE) and lie in [0,1]")
        reservoir_kind = 'specified_nonthermal'
    lead_in, lead_out = f[:, :, None, None]*gamma, (1-f)[:, :, None, None]*gamma
    lead_total = lead_R.sum(axis=0)
    shape = (len(energies), *H.shape)
    inside, outside, sigma = (np.zeros(shape, complex) for _ in range(3))
    shift = np.zeros_like(H)
    converged = False
    for iteration in range(1, max_iterations+1):
        GR = np.linalg.inv(energies[:, None, None]*np.eye(len(H))-H-shift-lead_total-sigma)
        GA = adjoint(GR)
        occupied = hermitian(GR@(lead_in.sum(axis=0)+inside)@GA)
        empty = hermitian(GR@(lead_out.sum(axis=0)+outside)@GA)
        new_in, new_out, new_R = phonon_self_energies(occupied, empty, M, bath, spacing)
        density = hermitian(np.trapezoid(occupied, energies, axis=0)/(2*np.pi))
        new_shift = hartree_shift(density, M, bath) if include_hartree else np.zeros_like(H)
        residual = max(float(np.max(abs(new_in-inside))), float(np.max(abs(new_out-outside))),
                       float(np.max(abs(new_R-sigma))), float(np.max(abs(new_shift-shift))))/model.t
        if residual < tolerance:
            converged = True
            break
        inside += mixing*(new_in-inside)
        outside += mixing*(new_out-outside)
        sigma += mixing*(new_R-sigma)
        shift += mixing*(new_shift-shift)
    if not converged:
        raise RuntimeError(f"SCBA did not converge: residual {residual:.3e} after {iteration} steps")
    if basis is not None:
        GR, occupied, empty, lead_R, lead_in, lead_out, sigma, inside, outside, shift = (
            basis@a@adjoint(basis) for a in [GR, occupied, empty, lead_R, lead_in, lead_out, sigma, inside, outside, shift])
    return Solution(model, bath, energies, tuple(chemical_potentials), temperature, reservoir_kind, GR, occupied, empty,
                    lead_R, lead_in, lead_out, sigma, inside, outside, shift, iteration, residual)


def elastic_transmission(model, energies):
    """Landauer charge transmission and physical-spin numerator; bath omitted."""
    H, _ = model.operators()
    sigma, gamma = model.embeddings(np.asarray(energies))
    GR = np.linalg.inv(np.asarray(energies)[:, None, None]*np.eye(len(H))-H-sigma.sum(axis=0))
    transmitted = gamma[1]@GR@gamma[0]@adjoint(GR)
    charge = np.trace(transmitted, axis1=-2, axis2=-1).real
    spin = np.trace(np.kron(np.eye(model.sites), SZ)@transmitted, axis1=-2, axis2=-1).real
    return charge, spin


def scalar_gauge_transmission(model, energies):
    """Independent scalar reduction of a nearest-neighbor SU(2) wire."""
    H, _ = model.operators()
    scalar = np.eye(model.sites)*2*model.t
    for j in range(model.sites-1):
        block = H[2*j:2*j+2, 2*j+2:2*j+4]
        hopping = np.sqrt(np.trace(block@block.conj().T).real/2)
        scalar[j, j+1] = scalar[j+1, j] = -hopping
    lead_R, gamma = model.embeddings(np.asarray(energies))
    sigma = lead_R.sum(axis=0)[:, ::2, ::2]
    GR = np.linalg.inv(np.asarray(energies)[:, None, None]*np.eye(model.sites)-scalar-sigma)
    return 2*gamma[0, :, 0, 0].real*gamma[1, :, -1, -1].real*abs(GR[:, -1, 0])**2


def dilute_scattering(model, energy, phonon_energy=.4, states=4, source=0, return_retarded=False):
    """Exact truncated one-electron, ground-phonon multichannel scattering.

    No Fermi sea is modeled here. This is an independent dilute/vacuum benchmark,
    not a replacement for the finite-occupation SCBA calculation.
    """
    if not isinstance(states, (int, np.integer)) or states < 2 or source not in (0, 1) or not np.isfinite([energy, phonon_energy]).all() or phonon_energy <= 0:
        raise ValueError("at least two phonon states, positive energy, and source 0/1 required")
    H, M = model.operators()
    dimension = len(H)
    b = np.diag(np.sqrt(np.arange(1, states)), 1)
    total = np.kron(np.eye(states), H)+np.kron(np.diag(np.arange(states)*phonon_energy), np.eye(dimension))+np.kron(b+b.T, M)
    sigma, gamma = model.embeddings(energy-np.arange(states)*phonon_energy)
    embedding = np.zeros_like(total)
    for n in range(states):
        block = slice(n*dimension, (n+1)*dimension)
        embedding[block, block] = sigma[:, n].sum(axis=0)
    GR = np.linalg.inv(energy*np.eye(len(total))-total-embedding)
    width_in = gamma[source, 0, 0 if source == 0 else -1, 0 if source == 0 else -1].real
    block_in = slice(0, 2) if source == 0 else slice(dimension-2, dimension)
    charge, spin, reflection, unitarity = 0., 0., 0., np.zeros((2, 2), complex)
    for port in (0, 1):
        for n in range(states):
            offset = n*dimension+(0 if port == 0 else dimension-2)
            block_out = slice(offset, offset+2)
            width = gamma[port, n, 0 if port == 0 else -1, 0 if port == 0 else -1].real
            amplitude = -1j*np.sqrt(width*width_in)*GR[block_out, block_in]
            if port == source and n == 0:
                amplitude += IDENTITY
            unitarity += amplitude.conj().T@amplitude
            flux = np.trace(amplitude@amplitude.conj().T).real/2
            if port != source:
                charge += flux
                spin += np.trace(SZ@amplitude@amplitude.conj().T).real/2
            else:
                reflection += flux
    result = dict(transmission=charge, reflection=reflection, transmitted_spin=spin,
                  transmitted_Pz=spin/charge if charge > 1e-15 else 0.,
                  unitarity_error=float(np.linalg.norm(unitarity-IDENTITY)))
    if return_retarded:
        result["ground_retarded"] = GR[:dimension, :dimension]
    return result
