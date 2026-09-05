# CISST1

Reproducible calculations for **Relaxation and Reciprocity in Chiral Spin
Transport**, by Mohamad Niknam and Louis-S. Bouchard (September 4, 2026).

Repository: https://github.com/lsbouchard/CISST1

The minimal Overleaf ZIP is separate: it contains only self-contained `main.tex`
and the seven figure PDFs actually used by the paper. Bibliography and tables
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
.venv/bin/python scripts/validate_outputs.py
.venv/bin/python scripts/test_models.py
```

The generators recreate eight figure PDFs, all CSV/JSON data, and two tables.
Do not edit generated numbers. Change model inputs in the scripts and rerun.
PDF timestamps are stripped. Repeatability is tested on the same pinned runtime;
other BLAS libraries can change last digits, so validation uses tolerances.

## Models

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

device = Device(cells=120)
Y, states = device.response(omega=1.0)
short_rates = device.decay_rates(load=0.0)
open_rates = device.decay_rates(load=float("inf"))
rate = slowest_pole(length=1.0, D=1.0, gamma=0.2, kappa=1.0)
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

## Scientific Assumptions

The device assumes a scalar spin mode, constant coefficients, Markovian bulk
relaxation, passive contacts, and charge response faster than the spin bandwidth.
The conversion vertex is an input, not derived from chirality or trajectory
curvature. Physical spin, orbital angular momentum, total electronic angular
momentum, mechanical rotation, and Kramers pseudospin are distinct.

The same-device cross coefficients satisfy `Y_cL = -Y_Lc` with the stated port
directions. Symmetric contacts and an even source give common spin conversion
while the differential response cancels. A passive resistive electrical load
adds a positive-semidefinite term to the decay operator. This prediction does
not cover circuit memory, finite-bias nonlinear response, or contact/material
changes with load.

The thickness data are **synthetic**, with seed 20260904 and 2% Gaussian noise
in log rate. No experimental amplitude or fitted experimental relaxation time
is claimed. Fit uncertainties are conditional on the model and known noise.

## Verification

`validate_outputs.py` recalculates every series and both tables and rejects
nonfinite residuals and incomplete samples. `test_models.py` also checks an
independent boundary-value solver, limiting cases, extreme contacts, corrupted
data, power balance, reciprocity, passivity, mesh convergence, load-dependent
eigenvalue ordering, and identifiability.
The GitHub Actions workflow runs the complete regeneration and checks on Python
3.12 and 3.14. Its result is separate from the author-side PDF build.

The author workspace additionally has the BibTeX manuscript and a `make pdf`
target. The minimal Overleaf archive needs pdfLaTeX and REVTeX4-2, not Python or
BibTeX. No redistribution license has yet been selected for the public code.
