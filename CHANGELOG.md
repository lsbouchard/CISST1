# Change log

## Contact-resolved response release - 2026-09-04

- Retitled the paper `Relaxation and Reciprocity in Chiral Spin Transport`.
- Replaced the difference-only spin port by independent inward contact currents.
  Derived a passive spatial three-port admittance and its Onsager identity.
- Added free-energy power balance, common/differential drive selection, exact
  Robin poles, synthetic thickness inference, and electrical-load feedback.
- Separated equilibrium linear response from finite drift and finite-bias
  differential response; corrected the Kubo sign and electrostatic-bias parity.
- Fixed the isolated-doublet susceptibility and stated virtual-transition
  corrections; restricted the Pauli Hamiltonian to weak magnetic fields.
- Clarified static-potential torque and conjugate-port sign conventions.
  Separated volume/voltage and Galerkin/eigenfunction notation.
- Stabilized contact exponentials and roots; added a weighted uniqueness
  argument, independent BVP checks, and nonfinite/corrupt-data regressions.
- Added the reciprocal device generator, complete output validation, synthetic
  inference metadata/table, ten numerical regression tests, and a GitHub Actions
  workflow for Python 3.12 and 3.14.
- Corrected Bloom et al.'s authors, updated Cho et al.'s publication details,
  and added Eckvahl, Latawiec, and Chiesa radical-pair work without equating
  those correlated dynamics with scalar transport relaxation.
- Credited earlier relaxation/diffusion/Edelstein models and narrowed novelty
  to the complete device, operational tests, and inference.
- Removed the kinematic diffusion map from the paper, retaining code/data;
  improved plot labels and added full-width operational figures.
- Preserved two-column REVTeX, UCLA affiliation, no acknowledgments, and the
  GitHub-only data statement. Kept minimal Overleaf and calculation ZIPs separate.

## Angular-momentum distinction release - 2026-08-07

- Shortened the title to `Relaxation-Limited Chiral Barnett Response`.
- Made the physical current-induced electronic-angular-momentum density the
  primary response and retained the Barnett frequency only as an optional
  response coordinate.
- Stated the local-quasiequilibrium condition required to infer spin from that
  auxiliary coordinate; otherwise a direct current-to-spin coefficient is
  required.
- Added Hahn, Tenn, and Augustine's zero-field Barnett analysis and used it to
  state the rotor-coupling-reaction-torque criterion for any physical Barnett
  interpretation.
- Recast $\Omega_{\chi}$ and $H_{\chi}^{\rm aux}$ as auxiliary response
  coordinates; the electric field or reservoir bias remains the physical drive.
- Added a typed three-sector angular-momentum balance in which spin-orbit and
  bath-transfer torques cancel internally and the circuit supplies the net
  source.
- Derived the vector Markov closure by explicit block-generator adiabatic
  elimination, including all block dimensions, its Schur complement, its
  stability condition, and the leading memory correction.
- Proved that target-density and volumetric-rate drives are non-identifiable at
  fixed `T1`, because only their combined source enters the bulk equation.
- Added finite-device mode poles and an explicitly normalized Green-Kubo
  derivation of the charge-spin Onsager-Casimir relation.
- Added a prepared-order initial-condition protocol, while stating that coherent
  oscillations require a resolved coupled model and are not a universal CISS
  prediction or a one-$T_1$ response.
- Distinguished collective rigid rotation of the nuclear framework or lattice,
  electronic orbital angular momentum, physical electron spin, total electronic
  angular momentum, internal angular-momentum channels, and Kramers pseudospin.
- Reserved $\Omega_{\rm mech}$ for mechanical rotation, $\Omega_{\rm frame}$
  for a coordinate-frame rate, $\Omega_{\rm traj}$ for the helical trajectory
  estimate, and $\Omega_{\chi}$ for the nonequilibrium generalized force.
- Replaced language suggesting that projection converts orbital angular
  momentum into spin. The rank-two projector now only represents distinct
  operators in one low-energy subspace; spin response is carried by an explicit
  physical-spin--electronic-angular-momentum cross-susceptibility.
- Replaced the informal full-space/projected-matrix identification by an
  isometry, an explicit rank-2 projector on the microscopic Hilbert space, and
  a typed `2 x 2` compression map; added a basis-invariance check for every
  physical response tensor.
- Exhibited the exact affine reparameterization freedom between target-density
  and torque-rate sources and derived dissipative stability directly from the
  positive-definite symmetric part of the `3 x 3` relaxation generator.
- Stated throughout that the trajectory-based Barnett identification is a
  hypothetical diagnostic introduced in this manuscript, not an earlier theory
  or publication by the authors.
- Clarified the full rotation generator, eliminated the last unqualified
  angular-velocity symbols, and split the Coriolis identity to fit one column.
- Generalized the contact solver to independent left and right transparencies
  and both drift directions, and added differential, boundary, and integrated
  conservation residuals.
- Restricted the physical interpretation of total-flux Robin boundaries at
  finite drift, while retaining and validating the general analytic boundary
  problem; redesigned the displayed drift-free contact comparison for clear
  two-column rendering.
- Regenerated and validated all numerical outputs, rebuilt the self-contained
  REVTeX4-2 source, and inspected every final figure and PDF page.
- Rebuilt a minimal Overleaf archive containing only `main.tex` and figure PDFs;
  scripts, CSV data, validation code, and this change log are in the separate
  code archive.

## Audited release - 2026-07-30

- Converted the manuscript to two-column REVTeX4-2 `reprint` format and made
  Figure 5 fill one column.
- Preserved the central conclusion while separating three distinct open-system
  protocols: target-density relaxation, volumetric injection, and boundary
  injection.
- Corrected and typed the Pauli spin-orbit, spin-rotation, Barnett-field,
  drift-diffusion, boundary-flux, and Onsager-Casimir equations.
- Defined the rank-two projector, every projected `2 x 2` operator, all `3 x 3`
  response maps, thermodynamic ports, dimensions, parities, and assumptions.
- Expanded Appendix C with explicit time-reversal, enantiomer, projection, and
  canonical-correlation derivations in continuous prose.
- Recast the novelty relative to orbital Edelstein and relaxon theories around
  protocol-resolved dynamics, finite-length response, frequency dependence,
  contacts, and direct-inverse reciprocity.
- Updated and live-checked the 2024-2026 literature; recent preprints remain
  identified as preprints, and no uncited data agreement is claimed.
- Added the analytic finite-contact calculation, its CSV output, residual
  checks, and absorbing/partial/reflecting contact figure.
- Regenerated every figure and numerical table using the requested constants;
  added `scripts/validate_outputs.py` to check every generated value.
- Removed vector-rendering seams from Figure 5 and manually inspected all six
  figure files and all 12 final manuscript pages.
- Removed the acknowledgments section and set the data-availability URL to
  `https://github.com/lsbouchard/CISST1`.
- Verified clean builds with no undefined citations, undefined references, or
  overfull boxes, and confirmed visually that every deferred float rendered.
