#!/usr/bin/env python3
"""Independent regression tests for contacts and reciprocal device dynamics."""

import unittest
from unittest.mock import patch
import numpy as np
from scipy.integrate import solve_bvp

from solve_contact_transparency import solve_profile, contact_observables, residuals
from reciprocal_device import Device, slowest_pole, asymmetric_pole, infer_parameters
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
        for d in [Device(cells=60), Device(kappa_left=.1, kappa_right=2., a=-.8, cells=60),
                  Device(kappa_left=.1, kappa_right=2., a=-.8, cells=60, geometry="local")]:
            c, K, B, direct = d.stiffness()
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

    def test_asymmetric_spectrum(self):
        for L in [.05, 1., 30.]:
            for kl, kr in [(1., .1), (.1, 1.), (0., 1.), (1., 1.)]:
                exact = asymmetric_pole(L, 1., .2, kl, kr)
                swapped = asymmetric_pole(L, 1., .2, kr, kl)
                self.assertAlmostEqual(exact, swapped)
                if kl == kr:
                    self.assertAlmostEqual(exact, slowest_pole(L, 1., .2, kl))
                errors = [abs(Device(length=L, D=1., gamma=.2, kappa_left=kl,
                                     kappa_right=kr, cells=n).decay_rates()[0]-exact)
                          for n in [60, 480]]
                self.assertLess(errors[-1]/exact, 3e-6)
                self.assertLess(errors[-1], errors[0]/10)

    def test_load_determinant_and_local_feedback(self):
        for geometry in ["lumped", "local"]:
            d = Device(cells=40, geometry=geometry, a=2., D=.1)
            c, K, B, direct = d.stiffness()
            baseline = d.decay_rates(count=40)
            Y0, _ = d.response(0.)
            for R in [.1, 1., 10., np.inf]:
                rates = d.decay_rates(R, count=40)
                self.assertTrue((rates >= baseline-1e-10).all())
                prediction = Y0[0, 0].real/d.G if np.isinf(R) else (1+R*Y0[0, 0].real)/(1+R*d.G)
                self.assertAlmostEqual(np.exp(np.log(rates/baseline).sum()), prediction, delta=1e-9)
            if geometry == "local":
                _, K0, _, _ = Device(cells=40, a=2., D=.1).stiffness()
                correction = K-K0
                np.testing.assert_allclose(correction@np.ones(40), 0., atol=1e-12)
                self.assertGreaterEqual(np.linalg.eigvalsh(correction).min(), -1e-12)

    def test_model_misspecification_and_capacity(self):
        from robustness import contact_fits, elimination_data
        _, fits = contact_fits()
        self.assertLess(fits[0]['D'], .4)
        self.assertLess(fits[0]['rms_log_mismatch'], .01)
        result = elimination_data()
        self.assertGreater(result['schur_rate']/result['exact_rates'][0], 100)
        self.assertLess(abs(result['corrected_rate']/result['exact_rates'][0]-1), 1e-4)

    def test_invalid_geometry_and_asymmetric_inputs(self):
        with self.assertRaises(ValueError):
            Device(geometry="unknown").response(0.)
        for contacts in [(np.nan, 1.), (-1., 1.), (np.inf, 1.)]:
            with self.assertRaises(ValueError):
                asymmetric_pole(1., 1., 1., *contacts)

    def test_weighted_kinetic_projection(self):
        weights = np.array([2., 2., 3., 3.])
        spin = np.array([1., -1., .5, -.5])
        velocity = np.array([3., -3., -1., 1.])
        W = np.diag(weights)
        capacity = spin@W@spin
        P = np.outer(spin, weights*spin)/capacity
        np.testing.assert_allclose(P@P, P, atol=1e-14)
        np.testing.assert_allclose(W@P, P.T@W, atol=1e-14)
        np.testing.assert_allclose(P@np.ones(4), 0., atol=1e-14)
        charge = -1.
        beta = charge*spin@W@velocity
        h, field = .3, .2
        current = charge*velocity@W@(spin*h)
        source = spin@W@(charge*velocity*field)
        self.assertAlmostEqual(current, beta*h)
        self.assertAlmostEqual(source, beta*field)
        self.assertAlmostEqual(current*field-h*source, 0.)

    def test_frequency_elimination_remainder(self):
        from robustness import elimination_data
        r = elimination_data()
        G = np.array(r['generator_per_second'])
        for omega in [.001, .01, .1]:
            exact = G[0, 0]-1j*omega-G[0, 1]*G[1, 0]/(G[1, 1]-1j*omega)
            reduced = r['schur_rate']-1j*omega*r['capacity']
            bound = omega**2*abs(G[0, 1]*G[1, 0])/G[1, 1]**3/(1-omega/G[1, 1])
            self.assertLessEqual(abs(exact-reduced), bound+2e-11)


class SpectralTests(unittest.TestCase):
    def test_modal_reconstruction_and_port_factorization(self):
        from spectral_tests import spectrum
        E = np.diag([1., -1., -1.])
        for geometry in ["lumped", "local"]:
            device = Device(cells=40, D=.1, a=1.2, kappa_left=1., kappa_right=.1, geometry=geometry)
            rates, overlaps = spectrum(device)
            _, _, _, direct = device.stiffness()
            for omega in [.01, 1., 100.]:
                reconstructed = np.diag(direct)+E@overlaps.T@(overlaps/(rates-1j*omega)[:, None])
                measured, _ = device.response(omega)
                np.testing.assert_allclose(reconstructed, measured, atol=1e-10)
            for d in overlaps:
                Z = E@np.outer(d, d)
                self.assertAlmostEqual(Z[0, 1]**2, -Z[0, 0]*Z[1, 1], delta=1e-10)
            weights = overlaps[:, 0]**2
            for t in [.01, 1., 10.]:
                for order in range(5):
                    self.assertGreaterEqual(np.sum(weights*rates**order*np.exp(-rates*t)), 0.)

    def test_load_slopes_and_interlacing(self):
        from spectral_tests import spectrum
        device = Device(cells=40, D=.1, a=1.2, kappa_left=1., kappa_right=.1)
        rates, overlaps = spectrum(device)
        # A one-sided second-order stencil only uses passive nonnegative loads.
        step = 1e-4
        r1 = device.decay_rates(step, count=6)
        r2 = device.decay_rates(2*step, count=6)
        slope = (-3*rates[:6]+4*r1-r2)/(2*step)
        np.testing.assert_allclose(slope, overlaps[:6, 0]**2, rtol=2e-5, atol=1e-8)
        for load in [.01, 1., np.inf]:
            loaded = device.decay_rates(load, count=device.cells)
            self.assertTrue((loaded >= rates-1e-10).all())
            self.assertTrue((loaded[:-1] <= rates[1:]+1e-10).all())

    def test_tail_bounds_random_spectra_and_units(self):
        from spectral_tests import tail_bracket
        rng = np.random.default_rng(20260905)
        for _ in range(80):
            rates = np.r_[.3, np.sort(rng.uniform(2., 30., 7))]
            weights = rng.uniform(.001, 2., 8)
            cutoff = .95*rates[1]
            tail = np.sum(weights[1:]/rates[1:])
            alpha = min(.4, .6*(cutoff-rates[0])/weights[0])
            lo, hi = tail_bracket(rates[0], weights[0], tail, cutoff, alpha)
            exact = np.linalg.eigvalsh(np.diag(rates)+alpha*np.outer(np.sqrt(weights), np.sqrt(weights)))[0]
            self.assertLessEqual(lo, exact+1e-12)
            self.assertLessEqual(exact, hi+1e-12)
            scaled = tail_bracket(rates[0]*1e9, weights[0]*1e21, tail*1e12, cutoff*1e9, alpha/1e12)
            np.testing.assert_allclose(np.asarray(scaled)/1e9, [lo, hi], rtol=1e-12)
        self.assertEqual(tail_bracket(1., 2., 0., 10., .5), (2., 2.))
        self.assertEqual(tail_bracket(1., 2., .1, 10., 0.), (1., 1.))
        for arguments in [(1., 2., .1, 1.1, 1.), (1., 2., -.1, 3., .1),
                          (1., 2., .1, np.nan, .1), (1., 0., .1, 3., .1)]:
            with self.assertRaises(ValueError):
                tail_bracket(*arguments)

    def test_dark_and_degenerate_modes(self):
        rates = np.array([.1, 1., 1., 4.])
        vertex = np.array([0., 2., 3., 1.])
        for alpha in [.1, 1., 10.]:
            loaded = np.linalg.eigvalsh(np.diag(rates)+alpha*np.outer(vertex, vertex))
            self.assertAlmostEqual(loaded[0], .1)
            self.assertAlmostEqual(loaded[1], 1.)
        # A degenerate level has a PSD sum of residues, not necessarily rank one.
        overlaps = np.array([[1., 2., 3.], [2., -1., 1.]])
        residue = overlaps.T@overlaps
        self.assertLess(residue[0, 1]**2, residue[0, 0]*residue[1, 1])

    def test_heterogeneous_local_elimination(self):
        from spectral_tests import series_correction
        widths = np.array([.1, .2, .3, .4])
        rho = np.array([1., 3., .4, 2.])
        beta = np.array([1., -.5, 2., -.8])
        G, vertex, K = series_correction(widths, rho, beta)
        h, voltage = np.array([.4, .3, -.2, .1]), .7
        current = G*voltage+vertex@h
        field = rho*(current-beta*h)
        self.assertAlmostEqual(widths@field, voltage)
        np.testing.assert_allclose(widths*beta*field, vertex*voltage-K@h, atol=1e-14)
        self.assertAlmostEqual(np.sum(widths*field**2/rho), G*voltage**2+h@K@h)
        self.assertGreaterEqual(np.linalg.eigvalsh(K).min(), -1e-12)
        np.testing.assert_allclose(K@(1/beta), 0., atol=1e-14)
        self.assertGreater(np.linalg.norm(K@np.ones(4)), .1)
        _, _, uniform = series_correction(np.ones(4)/4, np.ones(4), np.ones(4)*2)
        np.testing.assert_allclose(uniform, np.eye(4)-np.ones((4, 4))/4)
        for arrays in [([1.], [1., 2.], [1.]), ([1.], [0.], [1.]), ([1.], [1.], [np.nan])]:
            with self.assertRaises(ValueError):
                series_correction(*arrays)

    def test_multimode_kinetic_reduction(self):
        weights = np.diag([2., 3., 4.])
        modes = np.diag([1., 2., 3.])
        capacities = np.diag(modes.T@weights@modes)
        projector = sum(np.outer(modes[:, n], modes[:, n]@weights)/capacities[n] for n in [0, 1])
        np.testing.assert_allclose(projector@projector, projector)
        spin, velocity = np.array([1., -2., .5]), np.array([.3, .4, .6])
        h = np.array([.1, .2])
        psi = modes[:, :2]@h
        p = spin@weights@modes/capacities
        self.assertAlmostEqual(spin@weights@psi, np.sum(p[:2]*capacities[:2]*h))
        beta = -modes.T@weights@velocity
        self.assertAlmostEqual(-velocity@weights@psi, beta[:2]@h)
        fast_rates, fast_weights = np.array([10., 20.]), np.array([1., 3.])
        dc = np.sum(fast_weights/fast_rates)
        for omega in [.1, 1., 10.]:
            response = np.sum(fast_weights/(fast_rates-1j*omega))
            self.assertLessEqual(abs(response-dc), abs(omega)/min(fast_rates)*dc)


class HahnMechanismTests(unittest.TestCase):
    def test_rotating_hamiltonian_and_coupling(self):
        from hahn_transport import Helix, SY, SZ
        model = Helix()
        for time in [0., .3, 2., 5.]:
            U = model.unitary(time)
            a = model.zeta*model.omega
            lab = a*(-U@SY@U.conj().T+model.chi*model.pitch_ratio*SZ)/2
            np.testing.assert_allclose(U.conj().T@lab@U-model.rotation*SZ/2,
                                       model.hamiltonian, atol=1e-14)
            coupling_lab = U@model.coupling@U.conj().T
            np.testing.assert_allclose(U.conj().T@coupling_lab@U, model.coupling, atol=1e-14)
            np.testing.assert_allclose(U.conj().T@SZ@U, SZ, atol=1e-14)
            step = 1e-5
            plus, minus = model.unitary(time+step), model.unitary(time-step)
            derivative_phi = (plus@model.coupling@plus.conj().T
                              -minus@model.coupling@minus.conj().T)/(2*step*model.rotation)
            np.testing.assert_allclose(derivative_phi,
                                       -1j*(SZ@coupling_lab-coupling_lab@SZ)/2,
                                       atol=1e-10)
        mass, radius, pitch, omega = 2., .7, 1.3, -.8
        for chirality in [-1, 1]:
            axial_momentum = mass*chirality*pitch*omega
            orbital_lz = mass*radius**2*omega
            self.assertAlmostEqual(orbital_lz+chirality*pitch*axial_momentum,
                                   mass*(radius**2+pitch**2)*omega)

    def test_thermal_rates_and_gibbs_state(self):
        from hahn_transport import Helix, dissipator
        from scipy.linalg import expm
        for model in [Helix(), Helix(theta=2., omega=.2), Helix(omega=3., zeta=.3)]:
            for frequency in [.001, .1, 1., 10.]:
                self.assertAlmostEqual(model.spectrum(-frequency)/model.spectrum(frequency),
                                       np.exp(-frequency/(model.theta*model.omega_c)), delta=1e-13)
            down, up, _ = model.rates()
            gap = np.linalg.norm(model.vector)
            self.assertAlmostEqual(up/down, np.exp(-gap/(model.theta*model.omega_c)))
            rho = expm(-model.hamiltonian/(model.theta*model.omega_c))
            rho /= np.trace(rho)
            np.testing.assert_allclose(dissipator(rho, model.jumps()), 0., atol=1e-13)
            self.assertAlmostEqual(model.target(), -model.axis[2]*np.tanh(gap/(2*model.theta*model.omega_c)))

    def test_master_equation_and_positivity(self):
        from hahn_transport import Helix, IDENTITY, SZ, dissipator
        from scipy.integrate import solve_ivp
        model = Helix()
        down, up, _ = model.rates()
        total = down+up
        jumps = model.jumps()
        def rhs(t, vector):
            rho = vector.reshape((2, 2))
            H = model.hamiltonian
            return (-1j*(H@rho-rho@H)+dissipator(rho, jumps)).ravel()
        times = np.linspace(0, 4/total, 25)
        result = solve_ivp(rhs, (0, times[-1]), (IDENTITY/2).ravel(), t_eval=times,
                           rtol=1e-8, atol=1e-10, method="DOP853")
        self.assertTrue(result.success)
        for time, vector in zip(times, result.y.T):
            rho = vector.reshape((2, 2))
            self.assertAlmostEqual(np.trace(rho).real, 1., places=9)
            self.assertGreaterEqual(np.linalg.eigvalsh(rho).min(), -1e-9)
            self.assertAlmostEqual(np.trace(rho@SZ).real, float(model.polarization(time)), delta=2e-8)

    def test_frame_covariance_of_open_dynamics(self):
        from hahn_transport import Helix, SX, SZ, IDENTITY, dissipator
        from scipy.integrate import solve_ivp
        model = Helix()
        initial = (IDENTITY+.3*SX+.2*SZ)/2
        jumps = model.jumps()
        def rhs(t, vector, lab):
            rho = vector.reshape((2, 2))
            U = model.unitary(t)
            H = model.hamiltonian
            operators = jumps
            if lab:
                H = U@H@U.conj().T+model.rotation*SZ/2
                operators = [U@j@U.conj().T for j in jumps]
            return (-1j*(H@rho-rho@H)+dissipator(rho, operators)).ravel()
        times = np.linspace(0, 60., 31)
        results = [solve_ivp(lambda t, y: rhs(t, y, lab), (0, 60.), initial.ravel(),
                             t_eval=times, rtol=1e-9, atol=1e-11, method="DOP853") for lab in [False, True]]
        self.assertTrue(all(r.success for r in results))
        for i, t in enumerate(times):
            U = model.unitary(t)
            rotated = U@results[0].y[:, i].reshape((2, 2))@U.conj().T
            np.testing.assert_allclose(rotated, results[1].y[:, i].reshape((2, 2)), atol=2e-8)

    def test_zero_soc_chirality_and_flow(self):
        from hahn_transport import Helix
        times = np.array([0., 1., 100., 10000.])
        np.testing.assert_array_equal(Helix(zeta=0., g=0.).polarization(times), 0.)
        base = Helix().polarization(times)
        np.testing.assert_allclose(Helix(chi=-1).polarization(times), -base, atol=1e-13)
        np.testing.assert_allclose(Helix(flow=-1).polarization(times), -base, atol=1e-13)
        small = [float(Helix(zeta=.15*scale, g=.08*scale).polarization(100.)) for scale in [1e-3, 2e-3]]
        self.assertAlmostEqual(small[1]/small[0], 4., delta=.003)

    def test_lab_fixed_sidebands_and_zero_field_limit(self):
        from hahn_transport import Helix, SX, SP, SM
        for flow in [-1, 1]:
            model = Helix(zeta=0., flow=flow)
            down, up, _ = model.rates("lab_fixed")
            self.assertAlmostEqual(down, up, delta=1e-14)
            self.assertAlmostEqual(model.target("lab_fixed"), 0., delta=1e-14)
            t = .7
            U = model.unitary(t)
            np.testing.assert_allclose(U.conj().T@SX@U,
                SP*np.exp(1j*model.rotation*t)+SM*np.exp(-1j*model.rotation*t), atol=1e-14)
        self.assertGreater(abs(Helix().target()-Helix().target("lab_fixed")), .1)

    def test_residence_average_and_unit_scaling(self):
        from hahn_transport import Helix, flight_response
        from scipy.integrate import quad
        model = Helix()
        down, up, _ = model.rates()
        mean = 1000.
        integral = quad(lambda x: float(model.polarization(mean*x))*np.exp(-x), 0, np.inf)[0]
        expected = model.target()*mean*(down+up)/(1+mean*(down+up))
        self.assertAlmostEqual(integral, expected, places=10)
        self.assertLess(integral, float(model.polarization(mean)))
        scaled = Helix(omega=1e13, g=.08e13, omega_c=1e13)
        self.assertAlmostEqual(float(scaled.polarization(mean/1e13)), float(model.polarization(mean)), places=12)
        gamma = down+up
        for frequency in [0., gamma, 10*gamma, -gamma]:
            real = quad(lambda a: gamma*np.exp(-gamma*a)*np.cos(frequency*a), 0, mean)[0]
            imag = quad(lambda a: gamma*np.exp(-gamma*a)*np.sin(frequency*a), 0, mean)[0]
            self.assertAlmostEqual(flight_response(gamma, mean, frequency), real+1j*imag, places=12)
        self.assertAlmostEqual(flight_response(gamma, 100/gamma, gamma), 1/(1-1j), places=12)
        self.assertEqual(flight_response(0., mean, 0.), 0.)
        self.assertEqual(flight_response(gamma, 0., gamma), 0.)

    def test_invalid_inputs_and_weak_linewidths(self):
        from hahn_transport import Helix, example_data
        for kwargs in [dict(omega=0), dict(g=-1), dict(theta=0), dict(chi=0), dict(zeta=np.nan)]:
            with self.assertRaises(ValueError):
                Helix(**kwargs)
        with self.assertRaises(ValueError):
            Helix().rates("unspecified")
        with self.assertRaises(ValueError):
            Helix().polarization(-1.)
        _, _, rows, _ = example_data()
        for row in rows:
            self.assertLess(row["Gamma1_over_omega_c"], .01*row["omega_over_cutoff"])
            self.assertLess(abs(row["Pz_ten_turns"]), 1.)


if __name__ == "__main__":
    unittest.main(verbosity=2)
