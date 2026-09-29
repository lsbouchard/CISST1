# Finite Junction Transport Benchmark

September 29, 2026. This calculation tests whether the manuscript's real
Rashba-amplitude fluctuations can generate a longitudinal spin signal in a
finite, static helix connected to unpolarized electron reservoirs. It calculates
the outgoing flux rather than assigning it the bulk alignment target.

The reference calculation gives a small net drain spin signal, approximately
0.018%, under synthetic inputs. This is a conditional realization of
spin-selective transport, not a quantitative explanation of large experimental
CISS. Numerical convergence and conservation do not establish the molecular
coupling or the physical accuracy of a weak-coupling approximation.

## Model and Operator Spaces

The one-electron space is C^N tensor C^2, dimension D=2N, with basis ordered
as (site, physical spin up/down). The Pauli matrices act on physical electron
spin, not a Kramers pseudospin. Every electronic Hamiltonian, vertex,
embedding, Green matrix, and self-energy is D by D. Energy-sampled arrays
have shape (N_E,D,D); two-port arrays have shape (2,N_E,D,D).

Let u=z/b be the dimensionless axial coordinate, b the pitch divided by 2*pi,
chi=+/-1 the enantiomer sign, and phi_j=chi*u_j the actual azimuth. The ratio
r=R/b is dimensionless. The reference has N=9, angular step delta=pi/4,
r=0.7, signed mean SOC zeta=0.15, and fluctuation amplitude eta=0.08.
These are model inputs, not measured DNA parameters.

Define A_j=-sigma_phi(phi_j)+chi*r*sigma_z. The Hermitian lattice Hamiltonian
has on-site block 2t*1_2 and hopping block

```text
H[j,j+1] = -t*1_2 - i*t*zeta*delta*(A_j+A_{j+1})/4,
H[j+1,j] = H[j,j+1]^dagger.
```

This is a nearest-neighbor realization of the symmetrized continuum Rashba
interaction. Its finite lattice spacing is part of the benchmark model, not
a demonstrated continuum-limit error bound. In a continuum discretization,
t=hbar^2/(2*I_h*delta^2) and I_h=m_*(R^2+b^2), with units kg m^2. This is a
coordinate inertia, not the molecular rigid-body moment of inertia. The
relation zeta=2*alpha_R*m_*b/hbar^2 requires an independently established
effective mass and physical Rashba coefficient alpha_R [J m].

The bath vertex is exactly M=eta*dH/dzeta [energy]. The interaction is M tensor X,
where X=sum_a c_a(b_a+b_a^dagger) is dimensionless and Hermitian on the bosonic
bath space. The phonon operators act on the tensor product of harmonic-oscillator
Fock spaces. The retained phonon frequencies discretize a thermal Ohmic
spectrum; they do not truncate each oscillator's occupation in SCBA.

Both leads are spin-degenerate semi-infinite tight-binding chains, with on-site
energy 2t, hopping 2t, and spin-independent contact hopping t in the reference.
The retarded surface embedding is evaluated analytically. Its broadening is
positive semidefinite and vanishes outside the finite lead band. Reservoir
occupations are Fermi functions; the reference has mu_L=1.4t, mu_R=0.6t and
kB*T=0.15t. No externally rotating molecule, extra laboratory Barnett Zeeman
term, prescribed escape time, or imposed spin target appears in this H.

## Thermal Bath and Self Consistency

In reference units hbar=1, phonon energies are epsilon_a>0, occupations are
n_a=1/(exp(epsilon_a/T)-1), and w_a=c_a^2 are dimensionless positive weights.
The weights approximate the spectrum
S_X(+omega)=omega*exp(-omega/omega_c)/[omega_c^2*(1-exp(-omega/T))], with
omega_c=0.5t/hbar. The fine reference uses frequency spacing 0.025t/hbar and
maximum frequency 3t/hbar. Separate checks extend that cutoff to 4t/hbar.
The prescribed SI constants elsewhere in the repository are unchanged.

Write G^n=-i*G^< and G^p=i*G^>, both positive-semidefinite matrices with units
inverse energy. The full spatial and spin matrices enter the Fock self-energies:

```text
Sigma_ph^in(E)  = sum_a w_a M[(n_a+1)G^n(E+epsilon_a)
                             + n_a G^n(E-epsilon_a)]M,
Sigma_ph^out(E) = sum_a w_a M[(n_a+1)G^p(E-epsilon_a)
                             + n_a G^p(E+epsilon_a)]M.
Gamma_ph(E)    = Sigma_ph^in(E)+Sigma_ph^out(E).
```

Their retarded part is the causal principal-value transform of Gamma_ph,
minus i*Gamma_ph/2. The principal value is integrated exactly for a
piecewise-linear spectral interpolant with zero continuation beyond the
numerical window. The self-consistent Hartree displacement is
Sigma_H=-2*sum_a(w_a/epsilon_a)*Tr(M*rho)*M, with
rho=integral G^n(E)dE/(2*pi). The mean displacement and fluctuating vertex are
therefore not independent adjustable spin sources.

The Dyson and occupation equations are iterated together:

```text
G^R(E) = [E-H-Sigma_H-Sigma_leads^R(E)-Sigma_ph^R(E)]^(-1),
G^n(E) = G^R(E) Sigma_total^in(E) G^A(E),
G^p(E) = G^R(E) Sigma_total^out(E) G^A(E).
```

The unmixed fixed-point residual must be below 1e-11*t for published samples.
No diagonal self-energy approximation is used. FFT packing removes only
Hermitian redundancy; every spatial and spin off-diagonal element is retained.
At common chemical potential and temperature, the phonon collision terms obey
Sigma_ph^in(E)=f(E)*Gamma_ph(E). This establishes thermal balance and is tested
independently on arbitrary positive spectral matrices.

The conserving structure follows self-consistent electron-phonon NEGF, not an
added spin-alignment source. The thermal phonons act as an externally maintained
energy reservoir; no assertion of isolated electron-plus-phonon energy conservation
is made. [Charge-conservation derivation](https://iue.tuwien.ac.at/phd/pourfath/node61.html),
[coupled electron-phonon energy transport](https://arxiv.org/abs/0704.0723).

## Currents and Conservation

For port a, the particle current into the channel is the integral of
Tr[Sigma_a^in G^p-Sigma_a^out G^n]/(2*pi*hbar). Inserting physical sigma_z
gives a spin-weighted particle rate; multiplying that rate by hbar/2 gives
spin angular-momentum current. Conventional electrical current has electron
charge -e, which must not be confused with particle-flow direction.

Three spin observables are exported separately. The accumulated Pz is
Tr(sigma_z*rho)/Tr(rho). Gross extraction Pz is the spin imbalance of reservoir
extraction events. Net drain spin per particle subtracts both injection and
extraction. The net ratio is undefined at zero net particle current; CSVs leave
it blank rather than assigning a fictitious zero polarization. In general it
is not constrained by the bounds of a single outgoing spin density matrix.

Gross extraction is not source-resolved transmission: an electron extracted
at the right contact need not have been injected at the left contact. The
low-occupation drain control reduces this distinction, but does not establish
a general equivalence. A voltage or magnetoresistance detector would require
its own model.

Charge balance is checked between both leads and against the vanishing net
phonon number-collision integral. Energy balance is checked both from the
collision integral and independently from phonon emission minus absorption
events. The Hamiltonian and bath spin torques are retained in spin balance;
spin is not assumed conserved by SOC. The clamped helical structure and leads
can exchange orbital angular momentum with their supports. This benchmark
does not calculate a molecular rigid-body rotation or its reaction torque.

The spectral identity G^n+G^p=i(G^R-G^A), spectral sum rule, and one-particle
density eigenvalues are checked. The relative charge-balance residual is
normalized to the largest net or gross lead particle flux, so it remains
meaningful at equilibrium.
The lead-collision occupation formula assumes no independently occupied,
undamped bound state. A failed spectral sum rule therefore disqualifies a
parameter point rather than supplying a missing bound-state population.

## Independent Checks

For an elastic nearest-neighbor wire each hopping is a scalar magnitude times
a 2 by 2 SU(2) matrix. Site-dependent spin rotations remove those matrices on
an open chain. The spin-blind reservoir embeddings then give an unpolarized
elastic transmitted current. A separately constructed scalar wire reproduces
the charge transmission and Landauer current. This statement is specific to
this single-orbital-channel elastic model, not all multichannel CISS devices.

The Hamiltonian and fluctuation vertex are invariant under electron time
reversal. Opposite enantiomers are related by the physical-spin sigma_y
transformation. For symmetric contacts, spatial reflection combined with a
physical-spin rotation relates the two particle-flow directions. The numerical
reversal tests compare corresponding drain ports, not merely electrical
current sign conventions. Rotating every operator into the local spin frame
and transforming observables back leaves the physical solution unchanged.

An independent one-electron solver acts on C^(2N) tensor C^K, with K retained
phonon-number states, and diagonalizes the coupled-channel resolvent. It
explicitly sums elastic and inelastic reflection and transmission amplitudes.
Its full scattering probabilities sum to one, and phonon truncation is varied.
The leading perturbative Green correction agrees with exact scattering, with
the remaining error scaling as eta^4. A narrow, unpolarized, dilute incident
beam also compares SCBA charge and spin flux directly with this solver. The
beam test is explicitly nonthermal; production benchmark reservoirs are Fermi
baths. Agreement in this controlled limit does not certify SCBA at strong
coupling or high carrier density.

An independent change of energy unit rescales every Hamiltonian, lead, bath,
temperature, and chemical potential together. The density and polarization
remain unchanged, particle rates and spin torques scale linearly, and energy
currents scale quadratically, as required by the stated hbar=1 units.

## Resolution and Parameter Map

`data/junction/convergence.csv` records energy-grid and phonon-grid refinement,
phonon-cutoff extension, and electronic-window extension. The map records
the actual grid and bath resolution used for every point. Successive joint
energy/bath-grid differences must be below 3% for the net spin ratio and 0.2% for particle current,
with spectral sum-rule error below 0.003. These are numerical screening
criteria, not rigorous error bounds or uncertainties in a molecular prediction.
Cutoff and window extensions are tested separately for the reference, not for
every parameter point.

`feasibility.csv` varies geometrical length and eta. `variations.csv` varies
contacts, temperature, joint mean/fluctuating SOC scale, bias, and drain
occupation. `controls.csv` contains equilibrium, zero-SOC, elastic, enantiomer,
flow-reversal, and spin-frame controls. `dilute_scattering.csv` contains the
independent unitary scattering calculation. `parameters.json` records units,
the reference model, resolution comparisons, and the largest resolved map
signal. None of these files contains experimental data.

The stored elastic injectivity dwell time, evaluated at E=t, is only a reference residence scale.
It is not a directed active alignment time. Longer contact residence can
include backtracking, which changes the sampled rotation and its spin bias.
The steady-state self-energy width is not labeled T1: physical-spin
population relaxation requires a dynamical, vertex-consistent response
calculation and may have several modes.

## Resolved Results

For the fine reference, particle flux into the drain is 0.184766*t/hbar.
The spin observables are different quantities:

| Observable | Reference value |
| --- | ---: |
| Net drain spin per particle | 0.0177309% |
| Gross drain extraction polarization | 0.0111800% |
| Accumulated channel polarization | 0.00276051% |
| Heat delivered to phonons | 3.04559e-5*t^2/hbar |
| Maximum retarded phonon self-energy norm divided by t | 0.00145719 |

The medium/fine net-ratio difference is 0.00472%; extending the phonon cutoff
changes it by 0.00427%. Extending the electronic window changes it by less than
1e-8% relatively. These figures concern this reference observable, not a
uniform error bound on every density-matrix element.

All twelve length/noise map points pass the screening criteria. The largest
net drain ratio is 0.0708054% at one geometrical turn and eta=0.16. At eta=0.08,
half, one, and one-and-a-half turns give 0.00647378%, 0.0177300%, and 0.0169333%,
respectively. This finite-junction length dependence is not the
monotonic buildup of a prescribed unidirectional flight. The noise dependence
is approximately quadratic over this weak-coupling map; neither trend is a
universal molecular prediction.

The contact/temperature/SOC/bias controls give a largest tested net ratio of
0.0937223%, at twice the reference mean and fluctuating SOC. These bounded
parameter scans do not determine a global maximum. Reducing the drain chemical
potential to -3t gives 0.0845824%; its gross extraction and net ratio remain
separately exported rather than assumed identical.

## Relation to the Manuscript and Literature

The calculation preserves the intended mechanism: real helical Rashba coupling,
SOC fluctuations, orbital recoil, and nonequilibrium carrier transport in a
static structure. It supports a conditional relaxation-assisted spin signal.
It does not validate the large stationary trajectory target as the polarization
of a short molecular junction. The few-turn, weak-coupling benchmark is not
a substitute for the manuscript's active-turn restriction.
The dc signal alone does not identify a dominant Barnett-dressed alignment
mode or a unique T1; that identification requires the dynamical calculation.

Hahn, Tenn, and Augustine motivate transforming the physical rotating
interaction and retaining the coupling that transfers angular momentum.
Their nuclear-spin construction does not supply a molecular Rashba vertex,
phonon spectrum, effective mass, or electronic transit time.
[Chapter](https://link.springer.com/chapter/10.1007/3-540-32627-8_1).

Upadhyay and Levy's conserving weak-coupling calculation is an important
benchmark for approximation artifacts. Their transport signal uses a
magnetization-reversal protocol with a spin-polarized analyzer lead; our direct
spin flux uses two spin-degenerate leads. Their negligible signal in the
explored regime is not a numerical prediction for this vertex or observable,
and our small nonzero signal is not a contradiction of their result.
[Paper](https://arxiv.org/abs/2601.15063),
[journal DOI](https://doi.org/10.1063/5.0324706).

## Scientific Boundaries and Next Calculations

SCBA neglects higher-order vertex corrections and assumes phonons maintained
in a thermal state. Numerical self-energy size is a diagnostic, not a proof
of the approximation's accuracy. The common amplitude bath imposes spatially
correlated SOC noise. Real local phonons, disorder, hopping, and contacts can
produce different fluctuation vertices and scatter between helical sectors.
There is no explicit Coulomb charging term or bias-dependent electrostatic
potential. The stated Hartree shift is the phonon displacement, not an
electrostatic charging calculation. Chemical-potential bias changes reservoir
occupations while the bare channel Hamiltonian and contacts remain fixed.

No specific molecular junction has been identified or fitted here. Geometry,
mean SOC, the fluctuation vertex and phonon spectrum must be determined
independently for one system before the map becomes a material prediction.
The next discriminating calculation is a physical-spin dynamical response
with conserving vertex corrections, followed by comparison of collective and
local phonon couplings. Only then should a measured relaxation time be matched
to an alignment mode. Inverse CISS still requires the same reservoir-completed
model with explicitly conjugate spin and charge perturbations.

For an experimental comparison, independently constrain geometry, channel
energies, contact hopping, mean SOC, and the phonon modulation of SOC for one
junction family. Use those same inputs for its length, bias, temperature, and
enantiomer series, with an explicit spin-detector response. A fit with a separate
spin amplitude for every sample does not test this mechanism. A resolved
time-dependent spin response would additionally constrain relaxation; neither
the elastic dwell time nor the dc linewidth can substitute for that measurement.
Handedness and particle-flow reversals are necessary symmetry controls for the
specified symmetric-contact model, not a unique experimental signature of Barnett
alignment.

The current manuscript and minimal Overleaf bundle are not silently modified
by this research benchmark. Its quantitative result should be integrated in
the next scientific revision with these scope limits intact.

## Verification

`make validate` passes all 50 regression tests: 37 existing model tests and
13 new junction tests, including the deliberately shuffled figure-label check.
The data validator checks all 54 stored samples, including the 38 SCBA samples
and 16 independent scattering samples, and recomputes one SCBA map point from
scratch without the checkpoint. Every selected map/variation point has a
recorded resolution comparison. Both figure PDFs were rendered at 160 dpi and
manually inspected; the control-label order and a clipped title were corrected.

`make junction` regenerates the full calculation. `make junction-figures`
validates the stored data and redraws the figures without repeating the full
scan. The source hashes, numerical runtime, units, and solver settings are
recorded in `data/junction/parameters.json`. Distribution and GitHub verification
records are separate from the manuscript's unchanged compilation record.
