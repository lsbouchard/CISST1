#!/usr/bin/env python3
"""Validate generated numerical data against the manuscript equations."""

from pathlib import Path
import csv
import math


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

HBAR = 1.054571817e-34
KB = 1.380649e-23
GAMMA_E_ABS = 1.76085963023e11
N_TURN = 10.5


def read_rows(filename):
    with (DATA / filename).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(row, key):
    return float(row[key])


def assert_close(actual, expected, label, rtol=2e-12, atol=1e-14):
    if not math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol):
        raise AssertionError(
            f"{label}: actual={actual:.17g}, expected={expected:.17g}"
        )


def validate_constants():
    rows = read_rows("model_constants.csv")
    actual = {row["constant"]: float(row["value"]) for row in rows}
    expected = {
        "hbar": HBAR,
        "kB": KB,
        "abs_gamma_e": GAMMA_E_ABS,
        "N_turn": N_TURN,
    }
    if actual.keys() != expected.keys():
        raise AssertionError(f"constant names differ: {actual.keys()}")
    for key, value in expected.items():
        assert_close(actual[key], value, key, rtol=0.0, atol=0.0)


def validate_barnett():
    rows = read_rows("barnett_estimate.csv")
    if len(rows) != 601:
        raise AssertionError(f"expected 601 Barnett rows, found {len(rows)}")
    for index, row in enumerate(rows):
        rate = as_float(row, "R_hop_s_inv")
        omega = 2.0 * math.pi * rate / N_TURN
        assert_close(as_float(row, "Omega_rad_s_inv"), omega, f"Omega[{index}]")
        assert_close(
            as_float(row, "B_equiv_T"),
            omega / GAMMA_E_ABS,
            f"B_B[{index}]",
        )
        for temperature in (300.0, 100.0, 10.0):
            expected = math.tanh(HBAR * omega / (2.0 * KB * temperature))
            key = f"P_eq_T_{int(temperature)}"
            assert_close(as_float(row, key), expected, f"{key}[{index}]")

    selected = read_rows("barnett_selected_values.csv")
    if len(selected) != 6:
        raise AssertionError(
            f"expected 6 selected Barnett rows, found {len(selected)}"
        )
    for index, row in enumerate(selected):
        rate = as_float(row, "R_hop_s_inv")
        omega = 2.0 * math.pi * rate / N_TURN
        polarization = math.tanh(HBAR * omega / (2.0 * KB * 300.0))
        assert_close(as_float(row, "Omega_rad_s_inv"), omega, f"selected Omega[{index}]")
        assert_close(
            as_float(row, "B_equiv_T"),
            omega / GAMMA_E_ABS,
            f"selected B_B[{index}]",
        )
        assert_close(
            as_float(row, "P_eq_percent_300K"),
            100.0 * polarization,
            f"selected P_B[{index}]",
        )


def parse_tex_number(cell):
    """Parse the restricted scientific-number syntax used in the table."""

    token = cell.strip().strip("$")
    marker = r"\times10^{"
    if marker not in token:
        return float(token)
    mantissa, exponent = token.split(marker, maxsplit=1)
    return float(mantissa) * 10.0 ** int(exponent.rstrip("}"))


def validate_barnett_table():
    table_path = TABLES / "barnett_table.tex"
    rows = []
    for line in table_path.read_text(encoding="utf-8").splitlines():
        if "&" not in line or r"R_{\rm hop}" in line:
            continue
        cells = line.removesuffix(r"\\").split("&")
        rows.append([parse_tex_number(cell) for cell in cells])
    if len(rows) != 6 or any(len(row) != 4 for row in rows):
        raise AssertionError("Barnett table must contain six four-column data rows")

    for index, (rate, omega_shown, field_shown, percent_shown) in enumerate(rows):
        omega = 2.0 * math.pi * rate / N_TURN
        field = omega / GAMMA_E_ABS
        percent = 100.0 * math.tanh(HBAR * omega / (2.0 * KB * 300.0))
        assert_close(omega_shown, omega, f"table Omega[{index}]", rtol=6e-3)
        assert_close(field_shown, field, f"table B_B[{index}]", rtol=6e-3)
        assert_close(percent_shown, percent, f"table P_B[{index}]", rtol=6e-3)


def validate_relaxation():
    rows = read_rows("relaxation_response.csv")
    if len(rows) != 800:
        raise AssertionError(f"expected 800 relaxation rows, found {len(rows)}")
    for index, row in enumerate(rows):
        omega_t1 = as_float(row, "omega_T1")
        amplitude = 1.0 / math.sqrt(1.0 + omega_t1**2)
        lag = math.degrees(math.atan(omega_t1))
        assert_close(
            as_float(row, "normalized_amplitude"),
            amplitude,
            f"amplitude[{index}]",
        )
        assert_close(as_float(row, "lag_deg"), lag, f"lag[{index}]")


def validate_length_scaling():
    rows = read_rows("length_scaling.csv")
    if len(rows) != 700:
        raise AssertionError(f"expected 700 length rows, found {len(rows)}")
    for index, row in enumerate(rows):
        length = as_float(row, "normalized_length")
        drift = 1.0 - math.exp(-length)
        diffusion = 1.0 - 1.0 / math.cosh(length)
        flux = math.tanh(length)
        assert_close(
            as_float(row, "drift_exit_density_over_s_star"),
            drift,
            f"drift[{index}]",
        )
        assert_close(
            as_float(row, "diffusive_exit_density_over_s_star"),
            diffusion,
            f"diffusion[{index}]",
        )
        assert_close(
            as_float(row, "left_outward_flux_over_Ds_s_star_per_ell_s"),
            flux,
            f"flux[{index}]",
        )


def validate_diffusion_map():
    rows = read_rows("diffusion_landscape_grid.csv")
    if len(rows) != 200 * 220:
        raise AssertionError(f"expected 44000 diffusion rows, found {len(rows)}")
    for index, row in enumerate(rows):
        t1_ns = as_float(row, "T1_ns")
        diffusion = as_float(row, "D_s_m2_s")
        ell_nm = math.sqrt(diffusion * t1_ns * 1e-9) * 1e9
        assert_close(as_float(row, "ell_s_nm"), ell_nm, f"ell_s[{index}]")


def validate_contacts():
    validation = read_rows("contact_validation.csv")
    expected_validation_cases = {
        "absorbing",
        "partial",
        "reflecting",
        "asymmetric_positive_drift",
        "asymmetric_negative_drift",
        "reflecting_with_drift",
    }
    if {row["case"] for row in validation} != expected_validation_cases:
        raise AssertionError("contact validation cases are incomplete")
    for row in validation:
        if as_float(row, "max_abs_ode_residual") > 1e-10:
            raise AssertionError(f"ODE residual failed: {row}")
        if as_float(row, "max_abs_boundary_residual") > 1e-10:
            raise AssertionError(f"boundary residual failed: {row}")
        if as_float(row, "max_abs_balance_residual") > 1e-10:
            raise AssertionError(f"balance residual failed: {row}")

    rows = read_rows("contact_transparency.csv")
    if len(rows) != 3 * 220:
        raise AssertionError(f"expected 660 contact rows, found {len(rows)}")
    for index, row in enumerate(rows):
        length = as_float(row, "L_over_ell_s")
        average = as_float(row, "average_s_over_s_star")
        left_flux = as_float(row, "left_outward_flux_norm")
        right_flux = as_float(row, "right_outward_flux_norm")
        total_flux = as_float(row, "total_outward_flux_norm")
        assert_close(total_flux, left_flux + right_flux, f"flux sum[{index}]")
        assert_close(total_flux, length * (1.0 - average), f"balance[{index}]")
        if row["case"] in {"partial", "reflecting"}:
            tau_left = as_float(row, "tau_left")
            tau_right = as_float(row, "tau_right")
            assert_close(
                left_flux,
                tau_left * as_float(row, "left_boundary_s_over_s_star"),
                f"left Robin[{index}]",
            )
            assert_close(
                right_flux,
                tau_right * as_float(row, "right_boundary_s_over_s_star"),
                f"right Robin[{index}]",
            )
        elif row["case"] == "absorbing":
            assert_close(
                as_float(row, "left_boundary_s_over_s_star"),
                0.0,
                f"left absorbing[{index}]",
            )
            assert_close(
                as_float(row, "right_boundary_s_over_s_star"),
                0.0,
                f"right absorbing[{index}]",
            )
        if row["case"] == "reflecting":
            assert_close(average, 1.0, f"reflecting mean[{index}]")


def validate_figure_files():
    expected = {
        "fig_model_schematic.pdf",
        "fig_barnett_estimate.pdf",
        "fig_relaxation_response.pdf",
        "fig_length_scaling.pdf",
        "fig_contact_transparency.pdf",
        "fig_spin_diffusion_landscape.pdf",
    }
    actual = {path.name for path in FIGURES.glob("*.pdf")}
    if actual != expected:
        raise AssertionError(f"figure set differs: actual={sorted(actual)}")
    for filename in expected:
        if (FIGURES / filename).stat().st_size < 1_000:
            raise AssertionError(f"figure file is unexpectedly small: {filename}")


def main():
    validate_constants()
    validate_barnett()
    validate_barnett_table()
    validate_relaxation()
    validate_length_scaling()
    validate_diffusion_map()
    validate_contacts()
    validate_figure_files()
    print("Validated constants, table, figures, and all generated numerical data.")


if __name__ == "__main__":
    main()
