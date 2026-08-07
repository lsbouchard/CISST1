#!/usr/bin/env python3
"""Validate generated numerical data against the manuscript equations."""

from pathlib import Path
import csv
import math


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

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
    if {row["case"] for row in validation} != {
        "absorbing",
        "partial",
        "reflecting",
    }:
        raise AssertionError("contact validation cases are incomplete")
    for row in validation:
        if as_float(row, "max_abs_ode_residual") > 1e-10:
            raise AssertionError(f"ODE residual failed: {row}")
        if as_float(row, "max_abs_boundary_residual") > 1e-10:
            raise AssertionError(f"boundary residual failed: {row}")

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
        if row["case"] == "partial":
            tau = as_float(row, "tau_left")
            assert_close(
                left_flux,
                tau * as_float(row, "left_boundary_s_over_s_star"),
                f"left Robin[{index}]",
            )
            assert_close(
                right_flux,
                tau * as_float(row, "right_boundary_s_over_s_star"),
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
        elif row["case"] == "reflecting":
            assert_close(left_flux, 0.0, f"left reflecting[{index}]")
            assert_close(right_flux, 0.0, f"right reflecting[{index}]")


def main():
    validate_constants()
    validate_barnett()
    validate_relaxation()
    validate_length_scaling()
    validate_diffusion_map()
    validate_contacts()
    print("Validated constants and all generated numerical data.")


if __name__ == "__main__":
    main()
