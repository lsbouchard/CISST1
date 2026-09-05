#!/usr/bin/env python3
"""Independent regression tests for contacts and reciprocal device dynamics."""

import unittest
from unittest.mock import patch
import numpy as np
from scipy.integrate import solve_bvp

from solve_contact_transparency import solve_profile, contact_observables, residuals
from reciprocal_device import Device, slowest_pole, infer_parameters
import validate_outputs as validator


class ContactTests(unittest.TestCase):
    def test_independent_bvp(self):
        for L, tl, tr, pe in [(2., 1., 1., 0.), (3., .3, 2., 1.5), (3., 2., .3, -1.5), (2., 0., 0., 2.)]:
            x = np.linspace(0, L, 101)
            result = solve_bvp(lambda x, y: np.vstack((y[1], pe*y[1]+y[0]-1)),
                               lambda yl, yr: np.array([yl[1]-(pe+tl)*yl[0], yr[1]-(pe-tr)*yr[0]]),
                               x, np.zeros((2, len(x))), tol=1e-9, max_nodes=10000)
            self.assertTrue(result.success)
            profile = solve_profile(L, tl, tr, peclet=pe)
            np.testing.assert_allclose(profile.value(x), result.sol(x)[0], atol=1e-8)

    def test_overflow_regressions(self):
        for L, pe in [(1000., 0.), (8., 100.), (8., -100.)]:
            p = solve_profile(L, 1., 1., peclet=pe)
            self.assertTrue(np.isfinite(contact_observables(p)).all())
            self.assertLess(max(residuals(p, 1., 1., False)), 1e-9)
        self.assertAlmostEqual(solve_profile(1000., 1.).average(), .999, places=12)

    def test_invalid_inputs(self):
        for kwargs in [dict(length_over_ell=np.nan, tau_left=1), dict(length_over_ell=1, tau_left=np.nan),
                       dict(length_over_ell=1, tau_left=-1), dict(length_over_ell=1, tau_left=1, peclet=np.inf)]:
            with self.assertRaises(ValueError):
                solve_profile(**kwargs)

    def test_limits_and_contact_mode(self):
        p = solve_profile(2., 1.)
        mean, ul, ur, fl, fr = contact_observables(p)
        self.assertAlmostEqual(fl, (1-np.exp(-2))/2)
        self.assertAlmostEqual(fr, fl)
        self.assertAlmostEqual((fr-fl)/2, 0.)
        for L in [1e-4, .1, 1., 10., 1000.]:
            self.assertAlmostEqual(solve_profile(L, 0.).average(), 1.)
            mean = solve_profile(L, 0., absorbing=True).average()
            self.assertAlmostEqual(mean, 1-2*np.tanh(L/2)/L, delta=1e-10)

    def test_validator_rejects_nan_and_forgery(self):
        original = validator.read_rows
        for mode in ["nan", "forgery"]:
            def changed(filename):
                rows = [dict(r) for r in original(filename)]
                if mode == "nan" and filename == "contact_validation.csv":
                    rows[0]["max_abs_ode_residual"] = "nan"
                if mode == "forgery" and filename == "contact_transparency.csv":
                    for r in rows:
                        if r["case"] == "partial":
                            r["average_s_over_s_star"] = "1"
                            for key in ("left_boundary_s_over_s_star", "right_boundary_s_over_s_star",
                                        "left_outward_flux_norm", "right_outward_flux_norm", "total_outward_flux_norm"):
                                r[key] = "0"
                return rows
            with patch.object(validator, "read_rows", changed):
                with self.assertRaises(AssertionError):
                    validator.validate_contacts()


class ReciprocalTests(unittest.TestCase):
    def test_reciprocity_passivity_and_power(self):
        E = np.diag([1., -1., -1.])
        for d in [Device(cells=60), Device(kappa_left=.1, kappa_right=2., a=-.8, cells=60)]:
            c, diag, off, B, direct = d.matrices()
            K = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
            for w in [0., .1, 1., 10., 1000.]:
                Y, states = d.response(w)
                np.testing.assert_allclose(Y, E @ Y.T @ E, atol=1e-11)
                self.assertGreaterEqual(np.linalg.eigvalsh((Y+Y.conj().T)/2).min(), -1e-11)
                f = np.array([.7+.2j, -.4+.1j, .3-.5j])
                h = states @ f
                incoming_power = np.vdot(f, Y @ f).real
                loss = d.G * abs(f[0])**2 + np.vdot(h, K @ h).real
                for i, n in [(1, 0), (2, -1)]:
                    g = direct[i]
                    loss += g * (abs(f[i])**2 - 2 * (np.conj(h[n])*f[i]).real)
                self.assertAlmostEqual(incoming_power, loss, delta=1e-10)

    def test_common_differential_and_enantiomer(self):
        d = Device(cells=80)
        Y, _ = d.response(1.)
        flipped, _ = Device(a=-d.a, cells=80).response(1.)
        self.assertLess(abs(Y[0, 1]-Y[0, 2]), 1e-11)
        self.assertGreater(abs(Y[0, 1]+Y[0, 2]), .01)
        self.assertAlmostEqual(flipped[0, 1], -Y[0, 1])
        np.testing.assert_allclose(d.decay_rates(10), Device(a=-d.a, cells=80).decay_rates(10), atol=1e-10)

    def test_mesh_convergence_and_load_ordering(self):
        for kappa in [0., .3, 3., np.inf]:
            exact = slowest_pole(1., 1., 1., kappa)
            errors = [abs(Device(kappa_left=kappa, kappa_right=kappa, cells=n).decay_rates()[0]-exact)
                      for n in [40, 80, 160]]
            self.assertLess(errors[-1]/exact, 5e-5)
            if kappa != 0:
                self.assertLess(errors[-1], errors[0]/8)
        d = Device(cells=60)
        rates = np.array([d.decay_rates(R) for R in [0., .1, 1., 10., np.inf]])
        self.assertTrue((np.diff(rates, axis=0) >= -1e-10).all())

    def test_noiseless_identifiability(self):
        lengths = np.geomspace(.05, 30., 18)
        truth = [.8, .3, 1.2]
        rates = [slowest_pole(L, *truth) for L in lengths]
        fit, singular, _ = infer_parameters(lengths, rates)
        self.assertTrue(fit.success)
        np.testing.assert_allclose(np.exp(fit.x), truth, rtol=1e-7)
        self.assertGreater(singular[-1], 1.)
        scaled_fit, _, _ = infer_parameters(lengths * 1e-9, np.asarray(rates) * 1e9)
        np.testing.assert_allclose(np.exp(scaled_fit.x), np.asarray(truth) * [1e-9, 1e9, 1.], rtol=1e-7)

    def test_unidentifiable_fit_and_pole_extremes(self):
        with self.assertRaises(ValueError):
            infer_parameters(np.ones(6), np.ones(6))
        self.assertAlmostEqual(slowest_pole(1., 1., 1., 1e20), 1 + np.pi**2)
        self.assertAlmostEqual(slowest_pole(1., 1., 1., 1e-15), 1., delta=3e-15)
        with self.assertRaises(ValueError):
            Device(cells=10).decay_rates(count=11)


if __name__ == "__main__":
    unittest.main(verbosity=2)
