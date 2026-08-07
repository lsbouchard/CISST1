# Change log

## Angular-momentum distinction release - 2026-08-07

- Shortened the title to “Relaxation-Limited Chiral Barnett Response: An
  Open-System Theory of CISS.”
- Added Hahn, Tenn, and Augustine's zero-field Barnett analysis and used it to
  state the rotor-coupling-reaction-torque criterion for any physical Barnett
  interpretation.
- Recast $\Omega_{\chi}$ and $H_{\chi}^{\rm aux}$ as auxiliary response
  coordinates; the electric field or reservoir bias remains the physical drive.
- Added a typed three-sector angular-momentum balance in which spin-orbit and
  bath-transfer torques cancel internally and the circuit supplies the net
  source.
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
- Stated throughout that the trajectory-based Barnett identification is a
  hypothetical diagnostic introduced in this manuscript, not an earlier theory
  or publication by the authors.
- Clarified the full rotation generator, eliminated the last unqualified
  angular-velocity symbols, and split the Coriolis identity to fit one column.
- Regenerated and validated all numerical outputs, recompiled the self-contained
  REVTeX4-2 source, and manually inspected all six figures and all 13 PDF pages.
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
