#!/usr/bin/env python3
"""Validate generated numerical data against the manuscript equations."""

from pathlib import Path
import csv
import math
import json
import numpy as np


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
    if not (math.isfinite(actual) and math.isfinite(expected)) or not math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol):
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
    if len(validation) != len(expected_validation_cases) or {row["case"] for row in validation} != expected_validation_cases:
        raise AssertionError("contact validation cases are incomplete")
    for row in validation:
        for key in ("max_abs_ode_residual", "max_abs_boundary_residual", "max_abs_balance_residual"):
            value = as_float(row, key)
            if not math.isfinite(value) or value < 0:
                raise AssertionError(f"invalid nonnegative finite residual: {row}")
        if as_float(row, "max_abs_ode_residual") > 1e-10:
            raise AssertionError(f"ODE residual failed: {row}")
        if as_float(row, "max_abs_boundary_residual") > 1e-10:
            raise AssertionError(f"boundary residual failed: {row}")
        if as_float(row, "max_abs_balance_residual") > 1e-10:
            raise AssertionError(f"balance residual failed: {row}")

    rows = read_rows("contact_transparency.csv")
    if len(rows) != 3 * 220:
        raise AssertionError(f"expected 660 contact rows, found {len(rows)}")
    for case in ("absorbing", "partial", "reflecting"):
        lengths = sorted(as_float(r, "L_over_ell_s") for r in rows if r["case"] == case)
        if len(lengths) != 220:
            raise AssertionError(f"missing or duplicate contact samples: {case}")
        for i, value in enumerate(lengths):
            assert_close(value, .05 + (8. - .05) * i / 219, f"contact length grid {case}[{i}]")
    from solve_contact_transparency import solve_profile, residuals
    for row in validation:
        if row["case"] in {"absorbing", "partial", "reflecting"}:
            continue
        L = as_float(row, "length_over_ell")
        tl, tr, pe = (as_float(row, k) for k in ("tau_left", "tau_right", "peclet"))
        current = residuals(solve_profile(L, tl, tr, peclet=pe), tl, tr, False)
        if not all(math.isfinite(v) and 0 <= v <= 1e-10 for v in current):
            raise AssertionError(f"recomputed contact residual failed: {row['case']}")
    for index, row in enumerate(rows):
        length = as_float(row, "L_over_ell_s")
        average = as_float(row, "average_s_over_s_star")
        left_flux = as_float(row, "left_outward_flux_norm")
        right_flux = as_float(row, "right_outward_flux_norm")
        total_flux = as_float(row, "total_outward_flux_norm")
        if row["case"] not in {"partial", "reflecting", "absorbing"}:
            raise AssertionError(f"unknown contact case: {row['case']}")
        # An independent zero-drift hyperbolic solution checks actual profiles,
        # not just algebraic relations between values stored by the generator.
        t = math.tanh(length / 2)
        if row["case"] == "absorbing":
            expected_mean, expected_flux = 1 - 2 * t / length, t
        elif row["case"] == "reflecting":
            expected_mean, expected_flux = 1., 0.
        else:
            tau = as_float(row, "tau_left")
            assert_close(tau, as_float(row, "tau_right"), "symmetric plotted contacts")
            expected_flux = tau * t / (tau + t)
            expected_mean = 1 - 2 * expected_flux / length
        assert_close(as_float(row, "peclet"), 0., "plotted drift")
        assert_close(average, expected_mean, f"independent contact mean[{index}]", atol=2e-12)
        assert_close(left_flux, expected_flux, f"independent left flux[{index}]", atol=2e-12)
        assert_close(right_flux, expected_flux, f"independent right flux[{index}]", atol=2e-12)
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


def validate_device_outputs():
    from reciprocal_device import Device, slowest_pole, infer_parameters, node_rate
    device = Device()
    rows = read_rows("device_admittance.csv")
    if len(rows) != 160:
        raise AssertionError("expected 160 admittance samples")
    for row, omega in zip(rows, np.logspace(-3, 3, 160)):
        assert_close(as_float(row, "omega"), omega, "admittance frequency")
        Y, _ = device.response(omega)
        for key, value in {"Y_cL_real": Y[0, 1].real, "Y_cL_imag": Y[0, 1].imag,
                           "Y_Lc_real": Y[1, 0].real, "Y_Lc_imag": Y[1, 0].imag,
                           "common_inverse_abs": abs(Y[0, 1] + Y[0, 2]),
                           "differential_inverse_abs": abs(Y[0, 1] - Y[0, 2]),
                           "minimum_dissipation_eigenvalue": np.linalg.eigvalsh((Y + Y.conj().T)/2).min(),
                           "reciprocity_error": np.max(abs(Y - np.diag([1., -1., -1.]) @ Y.T @ np.diag([1., -1., -1.])))}.items():
            assert_close(as_float(row, key), value, key, atol=2e-11)
    rows = read_rows("device_load_decay.csv")
    if len(rows) != 90:
        raise AssertionError("expected 90 electrical loads")
    short = device.decay_rates()[0]
    for row, load in zip(rows, np.logspace(-3, 3, 90)):
        for key, value in {"RG": load*device.G, "spatial_rate": device.decay_rates(load)[0],
                           "short_rate": short, "node_rate": node_rate(device, load)}.items():
            assert_close(as_float(row, key), value, key, atol=2e-10)
    lengths = np.geomspace(.05, 30., 18)
    exact = np.array([slowest_pole(L, 1., .2, 1.) for L in lengths])
    observed = exact * np.exp(np.random.default_rng(20260904).normal(0., .02, len(lengths)))
    fit, singular, covariance = infer_parameters(lengths, observed)
    inferred = np.exp(fit.x)
    rows = read_rows("synthetic_thickness.csv")
    if len(rows) != 18:
        raise AssertionError("expected 18 synthetic thicknesses")
    for i, row in enumerate(rows):
        for key, value in {"length": lengths[i], "exact_rate": exact[i], "synthetic_rate": observed[i],
                           "fitted_rate": slowest_pole(lengths[i], *inferred), "log_standard_deviation": .02}.items():
            assert_close(as_float(row, key), value, key, rtol=1e-8)
    rows = read_rows("synthetic_parameter_recovery.csv")
    if [r["parameter"] for r in rows] != ["D", "Gamma", "kappa"]:
        raise AssertionError("incorrect inference parameter set")
    table = (TABLES / "synthetic_inference.tex").read_text()
    for i, row in enumerate(rows):
        for key, value in {"truth": [1., .2, 1.][i], "fitted": inferred[i],
                           "local_log_standard_error": np.sqrt(covariance[i, i])}.items():
            assert_close(as_float(row, key), value, key, rtol=1e-6)
        cells = " & ".join(f"{as_float(row, k):.3f}" for k in ["truth", "fitted", "local_log_standard_error"])
        if cells not in table:
            raise AssertionError("inference table does not match data")
    summary = json.loads((DATA / "device_parameters.json").read_text())
    if summary["device_parameters"] != device.__dict__ or summary["seed"] != 20260904 or summary["fit_success"] is not True:
        raise AssertionError("incorrect device metadata")
    np.testing.assert_allclose(summary["singular_values"], singular, rtol=1e-6)
    np.testing.assert_allclose(summary["fit_parameter_scales"], fit.parameter_scales, rtol=1e-12)
    np.testing.assert_allclose(summary["log_covariance"], covariance, rtol=1e-6)
    assert_close(summary["short_rate"], short, "short rate", atol=2e-10)
    assert_close(summary["open_rate"], device.decay_rates(np.inf)[0], "open rate", atol=2e-10)
    assert_close(summary["reduced_chi_squared"], sum(fit.fun**2)/15, "reduced chi squared", rtol=1e-8)


def validate_figure_files():
    expected = {
        "fig_model_schematic.pdf",
        "fig_barnett_estimate.pdf",
        "fig_relaxation_response.pdf",
        "fig_length_scaling.pdf",
        "fig_contact_transparency.pdf",
        "fig_spin_diffusion_landscape.pdf",
        "fig_device_inference.pdf",
        "fig_reciprocal_ports.pdf",
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
    validate_device_outputs()
    validate_figure_files()
    print("Validated constants, table, figures, and all generated numerical data.")


if __name__ == "__main__":
    main()
