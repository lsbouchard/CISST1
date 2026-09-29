# CISST1

Reproducible calculations for **A Barnett Mechanism for Chiral Spin Transport**,
by Mohamad Niknam and Louis-S. Bouchard (September 29, 2026).

Repository: https://github.com/lsbouchard/CISST1

The minimal Overleaf ZIP is separate: it contains only self-contained `main.tex`
and the thirteen figure PDFs actually used by the paper. Bibliography and tables
are embedded in the TeX file. The codebase contains no manuscript PDF or
third-party publications. An optional diffusion landscape remains reproducible
but is not a manuscript figure.

## Reproduce

Python 3.12 or later is required by the pinned dependencies:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/generate_figures.py
.venv/bin/python scripts/solve_contact_transparency.py
.venv/bin/python scripts/reciprocal_device.py
.venv/bin/python scripts/robustness.py
.venv/bin/python scripts/spectral_tests.py
.venv/bin/python scripts/hahn_transport.py
.venv/bin/python scripts/quantum_helix.py
.venv/bin/python scripts/validate_outputs.py
.venv/bin/python scripts/test_models.py
```

The generators recreate fourteen figure PDFs, all CSV/JSON data, and three tables.
Do not edit generated numbers. Change model inputs in the scripts and rerun.
PDF timestamps are stripped. Repeatability is tested on the same pinned runtime;
other BLAS libraries can change last digits, so validation uses tolerances.

## Models

The main mechanism is implemented in `hahn_transport.py`: a selected directed
helical trajectory samples a rotating Rashba interaction and Rashba-amplitude
fluctuations. Transforming both gives a stationary spin-plus-bath problem.
Weak-coupling thermal transition rates produce a dressed alignment target and
population relaxation time, followed by finite-transit outgoing spin polarization.
Hahn, Tenn, and Augustine (2006) motivate the rotating-interaction/relaxation
construction; the molecular coupling assumed here is not established by their
nuclear-spin calculation.

`Helix` accepts traversal frequency magnitude `omega`, signed mean SOC ratio `zeta`,
radius-to-pitch-parameter ratio `pitch_ratio`, fluctuating SOC rate `g`, bath
cutoff `omega_c`, temperature ratio `theta=kB*T/(hbar*omega_c)`, and signed
chirality/flow. Frequencies and `g` are in 1/s. The bath spectrum has units s.
The mean Hamiltonian and all jump operators are explicit 2-by-2 matrices;
`hamiltonian` returns H_rot/hbar in rate units. `rates()` returns the downward
and upward dressed-state rates and the normalized energy eigenvectors.
`polarization(time)` returns laboratory physical-spin Pz for unpolarized entry.
`target()` is undefined (NaN) when both rates vanish, whereas the finite-time
polarization correctly stays zero. Remove both mean and fluctuating SOC for the
zero-SOC limit. The API deliberately requires nonzero traversal frequency; it
does not extrapolate a selected-trajectory model into an equilibrium junction.

`bath="covariant"` uses a Rashba-amplitude coupling rotating with the sampled
helix. `bath="lab_fixed"` uses a distinct laboratory sigma_x coupling, including
both frequency-shifted sideband rates. These are different physical baths,
not different coordinate descriptions of the same bath. `jumps()` implements
the covariant model's complete secular dissipator, including dephasing.

`flight_response(gamma, transit, frequency)` returns the dimensionless complex
exit/target response for harmonic target modulation at fixed axis and population
rate, using exp(-i*frequency*t). Rates/frequencies are in 1/s and transit is in s.
Its transit numerator is retained; the single-pole result is only a long-flight
limit. Current modulation generally changes more than this target alone.

The three main figures and their CSV/JSON data cover alignment, the joint
mean/fluctuating SOC limit at finite transit time, and bath-dependent finite
transit versus stationary polarization. Parameters are synthetic reference
units, not experimental molecular estimates. The model gives a dissipative
spin polarizer, not unequal coherent elastic transmission eigenvalues.
Inverse CISS still requires matching orbital occupations and reservoir feedback;
it cannot be inferred by reversing this trajectory formula alone.

## Transport Extensions

`generate_figures.py` implements the stated constants, hypothetical
trajectory-substituted Barnett estimate, one-mode frequency response, and
closed-form drift/diffusion laws. A kinematic diffusion length does not establish
continuum validity at a molecular scale.

`solve_contact_transparency.py` solves `u'' - Pe*u' - u + 1 = 0` with independent
total-flux Robin contacts or two absorbing endpoints. Length is measured in
`sqrt(D*T1)`, `Pe = v*sqrt(T1/D)`, and `tau_i = kappa_i*sqrt(T1/D)`.
Reported contact fluxes point outward. Endpoint-anchored exponentials and stable
roots avoid overflow in the tested long-layer/high-drift regimes. Prescribed
drift is not automatically an equilibrium Onsager model.

`reciprocal_device.py` implements the equilibrium, zero-field, drift-free
three-port model. `Device` takes SI inputs: length [m], D [m^2/s], gamma [1/s],
susceptibility [J s^2/m^3], kappa [m/s], conversion a [C/m^3], and conductance per
area G [S/m^2]. Examples use synthetic dimensionless reference units, recorded
in `data/device_parameters.json`. The state has `cells` components.
`response(omega)` returns a complex 3-by-3 admittance and a `cells`-by-3 state
response, mapping `(V, h_L, h_R)` to `(j_c, I_L, I_R)`. Both spin currents point
inward, and `h_i = mu_spin_i/hbar` has frequency units.

```python
import sys
sys.path.insert(0, "scripts")
from reciprocal_device import Device, slowest_pole, infer_parameters

device = Device(cells=120, geometry="lumped")
Y, states = device.response(omega=1.0)
short_rates = device.decay_rates(load=0.0)
open_rates = device.decay_rates(load=float("inf"))
rate = slowest_pole(length=1.0, D=1.0, gamma=0.2, kappa=1.0)
local_device = Device(cells=120, geometry="local")
```

`slowest_pole` is the exact lowest symmetric Robin pole. `infer_parameters`
fits positive D, gamma, and kappa to a thickness series with known log-rate
uncertainty. It returns the SciPy fit, Jacobian singular values, and local
log-parameter covariance. Geometric-mean data scales make the fit independent
of the input units; parameter values are returned in those input units.
Dimensionless log parameters are bounded to [-12, 12], and the scaling factors
are recorded in `fit.parameter_scales`. Failed, boundary-limited, or numerically
rank-deficient fits are rejected. SI-unit invariance is regression-tested.
A full-rank Jacobian does not exclude poor conditioning, model error, parameter
drift across samples, or incorrect mode identification.

`geometry="lumped"` denotes reversible coupling of all spin elements to one
voltage port. `geometry="local"` denotes a uniform series conductor with
sigma = G*length and beta = a*length. The latter adds
`(beta**2/sigma)*(h - mean(h))` to spin loss, as required by local charge
conservation. Fast charge response alone does not make these geometries equal.
`stiffness()` returns capacity per cell, the full symmetric loss matrix, port
coupling matrix, and bare port conductances. Both geometries obey the same
passivity, reciprocity, and resistive-load determinant identities.

`asymmetric_pole` computes the lowest scalar Robin rate for unequal finite
contacts. It intentionally excludes local electrical feedback. It is validated
against a separate finite-volume discretization, including a reflecting end.

`robustness.py` generates contact-misspecification fits, restricted-length
information matrices, asymmetric port responses, comparisons of the two
electrical geometries, and the dynamic-capacity counterexample. It replaces
roundoff-only plots with physical contact asymmetry. The same-model noisy
recovery data remain available as a regression example, not a robustness claim.

## Quantum Completion and Applicability

`quantum_helix.py` constructs a static, symmetrized orbital-spin Rashba Hamiltonian
and its conserved screw generator Q=p_phi+Sz. `fourier_operators` returns
H/hbar [1/s], Q/hbar, and Sz/hbar on 2*(2*cutoff+1) Fourier-spin basis states.
`sector_hamiltonian` returns the exact 2-by-2 fixed-Q bulk fiber, including
orbital kinetic recoil. `sector_polarization` gives a pure +/-z preparation's
spin dynamics, with longitudinal relaxation, dephasing, and coherent precession.
`common_momentum_polarization` averages the two Q sectors of an unpolarized
beam with the same physical incident orbital momentum. It is a fixed-time
bulk spin calculation, not a finite-junction transmission solver.

The trajectory approximation needs I_h*abs(Omega)/hbar >> 1, not merely a large
equivalent field. The assumed free-electron-mass DNA-like geometry gives only
0.0668 at 1e13 assumed steps/s. The effective mass, physical channels, and
occupations must be established independently. Generated SI scale data use
CODATA 2022 m_e=9.1093837139e-31 kg as an illustrative mass, not a molecular fit.

`alignment_turns(fraction, rate_over_rotation)` computes the constant-flight
active-turn requirement. The reference rate needs 326 turns per T1 and 751
turns for 90% of its stationary target; ten turns give Pz=0.02914, not 0.96536.
Contact storage counts only while the same active helical generator persists.
The new figure and CSVs expose these two independent applicability conditions.
Generic spatial phonons can scatter between Q sectors and require a different
reservoir-resolved kinetic model. No inverse-CISS coefficient is inferred from
this conditional forward polarization.

## Scientific Assumptions

`spectral_tests.py` provides `spectrum(device)`, the modal port decomposition;
`tail_bracket(rate, weight, dc_tail, cutoff, alpha)`, a deterministic bracket
on the first electrically visible loaded rate; and `series_correction`, the
positive local electrical loss for piecewise-constant nonuniform coefficients.
The bracket takes rate/cutoff in 1/s, weight in S/m^2/s, residual dc conductance
in S/m^2, and alpha = R/(1+R*G) in ohm m^2. It raises an error when its sufficient
gap condition fails. The caller must independently establish the omitted-rate
cutoff; an absent fitted peak is not evidence of a gap. Measurement and model
uncertainty are not part of the returned interval.

The spectral example exports every modal weight, frequency sample, and load
bound to CSV, and records the inputs in `spectral_parameters.json`. The loaded
rates are calculated by direct diagonalization, separately from the bounds.
No measurement is fitted. The weighted local correction is independently tested
against charge conservation, the uneliminated field, and electrical power.

`Device` assumes a matched scalar spin mode and susceptibility, constant
coefficients, Markovian bulk relaxation, passive contacts, and a specified
electrical geometry with charge response faster than the spin bandwidth.
Fast-sector elimination generally changes dynamic capacity: a stable static
Schur complement is not by itself the observed inverse relaxation time.
The conversion vertex is an input, not derived from chirality or trajectory
curvature. Physical spin, orbital angular momentum, total electronic angular
momentum, mechanical rotation, and Kramers pseudospin are distinct.

The same-device cross coefficients satisfy `Y_cL = -Y_Lc` with the stated port
directions. Symmetric contacts and an even source give common spin conversion
while the differential response cancels. A passive resistive electrical load
adds a positive-semidefinite term to the decay operator. This prediction does
not cover circuit memory, finite-bias nonlinear response, or contact/material
changes with load.

For a finite model, the product of all loaded-to-short-circuit pole ratios is
`(1 + R*Ycc(0))/(1 + R*G)`. In a truly uniform single-node limit,
`lambda_open/lambda_short = Ycc(0)/G`. An apparently dominant measured pole does
not justify this last formula if other electrical spin modes carry dc weight.
The bare G must be independently calibrated or resolved in a suitable frequency
window. Neither identity uniquely identifies CISS.

More usefully, every simple short-circuit pole has initial load slope equal to
its electrical spectral numerator weight. These electrical self-response weights
are nonnegative; arbitrary spin-detector and cross-response weights need not be.
An isolated visible pole also gives a rank-one signed three-port numerator;
degeneracy replaces its pairwise factorization equality by an inequality.
Dark modes and mode crossings must be tracked explicitly. These statements
require a self-adjoint relaxation generator and a memoryless electrical load,
not merely passivity. The manuscript generalizes the construction to nonuniform
layers and to several kinetic distribution modes with physical-spin readout.

The thickness data are **synthetic**, with seed 20260904 and 2% Gaussian noise
in log rate. No experimental amplitude or fitted experimental relaxation time
is claimed. Fit uncertainties are conditional on the model and known noise.

## Verification

`validate_outputs.py` recalculates every series and all three tables and rejects
nonfinite residuals and incomplete samples. `test_models.py` also checks an
independent boundary-value solver, limiting cases, extreme contacts, corrupted
data, power balance, reciprocity, passivity, mesh convergence, load-dependent
eigenvalue ordering, and identifiability.
Tests additionally cover unequal-contact root convergence, local feedback,
both-geometry power balance, determinant identities, and misleading fits under
contact-model misspecification. Spectral tests include reconstruction,
residue factorization, load-slope/interlacing checks, 80 randomized tail-bound
cases, SI rescaling, dark/degenerate modes, nonuniform local power balance, and
multimode kinetic projection with a controlled fast conductivity background.
The full suite now has 37 tests. Mechanism tests check the exact spin-frame
transformation, thermal KMS balance and stationary Gibbs state, direct master
equation integration and positivity, laboratory/rotating open-system covariance,
zero SOC and reversal symmetries, lab-fixed sideband cancellation, residence-time
averaging, finite-flight frequency quadrature, screw-coordinate bookkeeping,
physical-unit rescaling, invalid inputs, and weak rate/gap ratios.
Additional independent quantum tests check Hermiticity, screw conservation,
time-reversal symmetry and zero equilibrium spin, exact Fourier-block matching,
the fluctuation vertex, matrix-exponential spin dynamics, common-momentum
preparation, and the finite-turn constraint. Signed Rashba coefficients and
overflow-safe large finite electrical loads have regression coverage.
The GitHub Actions workflow runs the complete regeneration and checks on Python
3.12 and 3.14. Its result is separate from the author-side PDF build.

The author workspace additionally has the BibTeX manuscript and a `make pdf`
target. The minimal Overleaf archive needs pdfLaTeX and REVTeX4-2, not Python or
BibTeX. No redistribution license has yet been selected for the public code.
