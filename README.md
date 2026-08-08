# CISST1 calculation code

This archive contains the scripts and source CSV data for the figures and
numerical checks in *Relaxation-limited chiral Barnett response*. The minimal
Overleaf archive is distributed separately and contains only its self-contained
`main.tex` and six figure PDFs.

## Contents

- `scripts/generate_figures.py`: model schematic, Barnett estimates,
  relaxation response, length scaling, diffusion map, CSV files, and table
- `scripts/solve_contact_transparency.py`: analytic finite-contact calculation
- `scripts/validate_outputs.py`: independent checks of generated numbers,
  table rows, and the complete figure set
- `data/`: generated source data and analytic validation residuals
- `requirements.txt`: pinned Python dependencies
- `CHANGELOG.md`: manuscript and calculation changes

## Reproduce

Python 3.11 or later is required. From the archive root, run:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/generate_figures.py
.venv/bin/python scripts/solve_contact_transparency.py
.venv/bin/python scripts/validate_outputs.py
```

The first two scripts regenerate `figures/*.pdf`, `data/*.csv`, and
`tables/barnett_table.tex`. The validator checks the requested physical
constants, every displayed Barnett-table value, the one-mode response,
finite-length laws, the diffusion grid, Robin boundary conditions,
differential-equation residuals, integrated angular-momentum balance, and the
presence and nonzero size of all six figure PDFs. A successful run ends with
`Validated constants, table, figures, and all generated numerical data.`

## Scientific Scope

The manuscript distinguishes mechanical rotation, electronic orbital angular
momentum, physical spin, total electronic angular momentum, other internal
angular-momentum reservoirs, and Kramers pseudospin. The current-induced
electronic response is the physical source; the symbol `Omega_chi` is only an
auxiliary conjugate coordinate.

The open-system dynamics separate target-density relaxation, volumetric rate
injection, boundary injection, and conversion of prepared internal order. A
resolved balance identifies external source and sink torques and makes
spin-orbit coupling an internal transfer. Coherent prepared-order conversion
requires a coupled kinetic model and is not represented by the one-`T1`
relaxation equation.

The public calculation repository is
`https://github.com/lsbouchard/CISST1`.
