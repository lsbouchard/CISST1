#!/usr/bin/env python3
"""Independent operator, reservoir, conservation and scattering checks."""
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
import numpy as np
from scipy.special import expit
from junction_transport import (Junction, ThermalBath, solve, energy_grid, adjoint,
    shifted, principal_value, phonon_self_energies, elastic_transmission,
    scalar_gauge_transmission, dilute_scattering)
from hahn_transport import SY, SZ


class JunctionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Junction(sites=5, contact_left=1.4, contact_right=1.4)
        cls.bath = ThermalBath(temperature=.25)
        cls.grid = energy_grid(spacing=.05)
        cls.solution = solve(cls.model, cls.bath, cls.grid, tolerance=1e-11)

    def test_operators_and_signed_chirality(self):
        model = self.model
        H, M = model.operators()
        T = np.kron(np.eye(model.sites), 1j*SY)
        np.testing.assert_allclose(H, adjoint(H), atol=1e-14)
        np.testing.assert_allclose(M, adjoint(M), atol=1e-14)
        np.testing.assert_allclose(T@H.conj()@adjoint(T), H, atol=1e-14)
        np.testing.assert_allclose(T@M.conj()@adjoint(T), M, atol=1e-14)
        delta = 1e-4
        H1, _ = replace(model, zeta=model.zeta+delta).operators()
        np.testing.assert_allclose((H1-H)/delta, M/model.eta, atol=2e-12)
        Y = np.kron(np.eye(model.sites), SY)
        reverse_H, reverse_M = replace(model, chi=-1).operators()
        np.testing.assert_allclose(Y@H@Y, reverse_H, atol=1e-14)
        np.testing.assert_allclose(Y@M@Y, reverse_M, atol=1e-14)

    def test_elastic_scalar_gauge_and_landauer(self):
        grid = self.grid
        for zeta in [-.4, 0., .4]:
            model = replace(self.model, zeta=zeta, eta=0.)
            charge, spin = elastic_transmission(model, grid)
            np.testing.assert_allclose(charge, scalar_gauge_transmission(model, grid), atol=2e-13)
            np.testing.assert_allclose(spin, 0., atol=2e-13)
        model = replace(self.model, eta=0.)
        result = solve(model, self.bath, grid)
        charge, _ = elastic_transmission(model, grid)
        fL, fR = [expit((mu-grid)/self.bath.temperature) for mu in (1.4, .6)]
        expected = np.trapezoid(charge*(fL-fR), grid)/(2*np.pi)
        self.assertAlmostEqual(result.observables()['particle_into_left'], expected, delta=1e-13)

    def test_thermal_kms_and_fft_convolution(self):
        rng = np.random.default_rng(9341)
        grid = energy_grid(-2., 3., .05)
        bath = ThermalBath((.2, .4), (.3, .7), .25)
        B = rng.normal(size=(len(grid), 4, 4))+1j*rng.normal(size=(len(grid), 4, 4))
        spectral = B@adjoint(B)
        f = expit((.7-grid)/bath.temperature)[:, None, None]
        occupied, empty = f*spectral, (1-f)*spectral
        V = np.diag([1., -1., .3, -.2])
        inside, outside, sigma = phonon_self_energies(occupied, empty, V, bath, .05)
        np.testing.assert_allclose(inside, f*(inside+outside), atol=5e-13)
        reference_in, reference_out = np.zeros_like(occupied), np.zeros_like(empty)
        for energy, weight, n in zip(bath.energies, bath.weights, bath.occupations):
            offset = int(round(energy/.05))
            reference_in += weight*((n+1)*shifted(occupied, offset)+n*shifted(occupied, -offset))
            reference_out += weight*((n+1)*shifted(empty, -offset)+n*shifted(empty, offset))
        np.testing.assert_allclose(inside, V@reference_in@V, atol=1e-13)
        np.testing.assert_allclose(outside, V@reference_out@V, atol=1e-13)
        np.testing.assert_allclose(1j*(sigma-adjoint(sigma)), inside+outside, atol=1e-13)

    def test_principal_value_against_causal_lorentzian(self):
        errors = []
        for spacing in [.1, .05]:
            grid = energy_grid(-50., 50., spacing)
            width = .7
            broadening = (2*width/(grid*grid+width*width))[:, None, None]
            transformed = principal_value(broadening)[:, 0, 0].real
            selected = abs(grid) < 2
            error = max(abs(transformed[selected]-grid[selected]/(grid[selected]**2+width**2)))
            errors.append(error)
        self.assertLess(errors[1], .002)
        self.assertGreater(errors[0]/errors[1], 3.)

    def test_charge_energy_spin_balance_and_positive_state(self):
        obs = self.solution.observables()
        self.assertLess(obs['charge_balance_relative'], 1e-9)
        self.assertLess(abs(obs['phonon_number_collision']), 1e-10)
        self.assertLess(abs(obs['energy_balance']), 1e-12)
        self.assertLess(abs(obs['event_energy_balance']), 2e-10)
        self.assertLess(abs(obs['spin_balance']), 1e-12)
        self.assertLess(obs['occupied_empty_identity_error'], 1e-11)
        self.assertGreaterEqual(obs['density_min_eigenvalue'], -1e-10)
        self.assertLessEqual(obs['density_max_eigenvalue'], 1.001)
        self.assertLess(obs['spectral_sumrule_error'], .003)
        self.assertLess(abs(obs['extraction_Pz_right']), 1.)
        self.assertGreater(abs(obs['net_spin_per_particle_right']), 1e-7)

    def test_energy_unit_rescaling(self):
        scale = 10.
        model = replace(self.model, t=scale*self.model.t,
                        lead_hopping=scale*self.model.lead_hopping,
                        contact_left=scale*self.model.contact_left,
                        contact_right=scale*self.model.contact_right)
        bath = replace(self.bath, energies=tuple(scale*e for e in self.bath.energies),
                       temperature=scale*self.bath.temperature)
        scaled = solve(model, bath, scale*self.grid,
                       chemical_potentials=(1.4*scale, .6*scale), tolerance=1e-11)
        np.testing.assert_allclose(scale*scaled.retarded, self.solution.retarded, atol=2e-11)
        np.testing.assert_allclose(scaled.density, self.solution.density, atol=2e-12)
        base, other = self.solution.observables(), scaled.observables()
        for key in ['particle_into_right', 'spin_sigma_into_right',
                    'hamiltonian_spin_torque', 'bath_spin_torque']:
            np.testing.assert_allclose(other[key], scale*base[key], atol=2e-12, rtol=1e-8)
        for key in ['heat_into_phonons', 'heat_from_events', 'energy_into_right']:
            np.testing.assert_allclose(other[key], scale**2*base[key], atol=2e-11, rtol=1e-8)
        self.assertAlmostEqual(other['net_spin_per_particle_right'],
                               base['net_spin_per_particle_right'], delta=2e-11)

    def test_equilibrium_and_zero_soc(self):
        equilibrium = solve(self.model, self.bath, self.grid, chemical_potentials=(1., 1.), tolerance=1e-11)
        obs = equilibrium.observables()
        for key in ['particle_into_left', 'particle_into_right', 'spin_sigma_into_left',
                    'spin_sigma_into_right', 'accumulated_Pz', 'heat_into_phonons']:
            self.assertLess(abs(obs[key]), 1e-10, key)
        self.assertIsNone(obs['net_spin_per_particle_right'])
        spin_free = solve(replace(self.model, zeta=0., eta=0.), self.bath, self.grid)
        self.assertLess(abs(spin_free.observables()['spin_sigma_into_right']), 1e-13)

    def test_enantiomer_and_actual_particle_flow_reversal(self):
        base = self.solution.observables()
        other = solve(replace(self.model, chi=-1), self.bath, self.grid, tolerance=1e-11).observables()
        self.assertAlmostEqual(other['particle_into_right'], base['particle_into_right'], delta=1e-12)
        self.assertAlmostEqual(other['spin_sigma_into_right'], -base['spin_sigma_into_right'], delta=1e-12)
        reverse = solve(self.model, self.bath, self.grid, chemical_potentials=(.6, 1.4), tolerance=1e-11).observables()
        self.assertAlmostEqual(reverse['particle_into_left'], base['particle_into_right'], delta=1e-12)
        self.assertAlmostEqual(reverse['extraction_Pz_left'], -base['extraction_Pz_right'], delta=1e-11)
        self.assertAlmostEqual(reverse['net_spin_per_particle_left'], -base['net_spin_per_particle_right'], delta=1e-11)

    def test_rotating_spin_basis_is_only_a_representation_change(self):
        phi = self.model.chi*(np.arange(self.model.sites)-(self.model.sites-1)/2)*self.model.step
        U = np.diag(np.column_stack([np.exp(-.5j*phi), np.exp(.5j*phi)]).ravel())
        rotated = solve(self.model, self.bath, self.grid, tolerance=1e-11, basis=U)
        np.testing.assert_allclose(rotated.retarded, self.solution.retarded, atol=2e-11)
        np.testing.assert_allclose(rotated.density, self.solution.density, atol=2e-12)
        self.assertAlmostEqual(rotated.observables()['spin_sigma_into_right'],
                               self.solution.observables()['spin_sigma_into_right'], delta=2e-12)

    def test_dilute_exact_scattering_and_born_limit(self):
        model = replace(self.model, sites=3, eta=.02)
        results = [dilute_scattering(model, 1.2, states=n) for n in [3, 5]]
        for result in results:
            self.assertLess(result['unitarity_error'], 1e-13)
            self.assertAlmostEqual(result['transmission']+result['reflection'], 1., delta=1e-13)
        self.assertAlmostEqual(results[0]['transmitted_Pz'], results[1]['transmitted_Pz'], delta=1e-8)
        errors = []
        for eta in [.01, .005]:
            weak = replace(model, eta=eta)
            H, M = weak.operators()
            sigma, _ = weak.embeddings(np.array([1.2, .8]))
            G0 = np.linalg.inv(1.2*np.eye(len(H))-H-sigma[:, 0].sum(axis=0))
            G1 = np.linalg.inv(.8*np.eye(len(H))-H-sigma[:, 1].sum(axis=0))
            perturbative = G0+G0@M@G1@M@G0
            exact = dilute_scattering(weak, 1.2, return_retarded=True)['ground_retarded']
            errors.append(np.linalg.norm(exact-perturbative))
        self.assertTrue(15 < errors[0]/errors[1] < 17)

    def test_scba_against_independent_dilute_scattering(self):
        grid = energy_grid(spacing=.01)
        occupation = 1e-6*np.exp(-.5*((grid-1.2)/.1)**2)
        model = replace(self.model, sites=3, eta=.02)
        result = solve(model, ThermalBath(temperature=0.), grid, tolerance=1e-12,
                       lead_occupations=np.array([occupation, np.zeros_like(occupation)]))
        charge = np.zeros_like(grid)
        spin = np.zeros_like(grid)
        for i in np.flatnonzero(occupation > 1e-14):
            exact = dilute_scattering(model, grid[i], states=5)
            charge[i], spin[i] = exact['transmission'], exact['transmitted_spin']
        current = 2*np.trapezoid(occupation*charge, grid)/(2*np.pi)
        spin_current = 2*np.trapezoid(occupation*spin, grid)/(2*np.pi)
        obs = result.observables()
        self.assertAlmostEqual(-obs['particle_into_right']/current, 1., delta=.001)
        self.assertGreater(abs(spin_current), 1e-14)
        self.assertAlmostEqual(-obs['spin_sigma_into_right']/spin_current, 1., delta=.02)
        self.assertEqual(result.reservoir_kind, 'specified_nonthermal')

    def test_input_errors_and_unconverged_results_rejected(self):
        for kwargs in [dict(sites=2), dict(contact_left=0.), dict(eta=-1.), dict(chi=0)]:
            with self.assertRaises(ValueError):
                Junction(**kwargs)
        for kwargs in [dict(energies=-1.+.03*np.arange(201)), dict(basis=np.zeros((10, 10)))]:
            with self.assertRaises(ValueError):
                solve(self.model, self.bath, **kwargs)
        with self.assertRaises(RuntimeError):
            solve(self.model, self.bath, self.grid, max_iterations=1)

    def test_control_plot_names_follow_bar_values(self):
        import benchmark_junction as benchmark
        from matplotlib.figure import Figure
        map_rows = [dict(eta=eta, geometrical_turns=turns, net_spin_per_particle_right=.001)
                    for eta in [.04, .08, .16] for turns in [.5, 1., 1.5]]
        variations = [dict(parameter=name, value=value, net_spin_per_particle_right=.001,
                           extraction_Pz_right=.0005)
                      for name, values in [('contact', [.7, 1., 1.4]),
                                           ('temperature', [.075, .15, .3]),
                                           ('SOC_scale', [0., 1., 2.])]
                      for value in values]
        controls = [dict(check=name, net_spin_per_particle_right=value)
                    for name, value in [('elastic_SOC_only', 0.), ('reference', .001),
                                        ('all_SOC_zero', 0.), ('opposite_enantiomer', -.001)]]
        figures = []
        def capture(figure, *args, **kwargs):
            figures.append(figure)
        with TemporaryDirectory() as directory, patch.object(benchmark, 'FIGURES', Path(directory)), \
                patch.object(Figure, 'savefig', autospec=True, side_effect=capture):
            benchmark.write_figures(map_rows, variations, controls)
        axis = figures[1].axes[0]
        np.testing.assert_allclose([bar.get_height() for bar in axis.patches], [.1, -.1, 0., 0.])
        self.assertEqual([tick.get_text() for tick in axis.get_xticklabels()],
                         ['Reference', r'$\chi=-1$', 'Zero SOC', 'Elastic'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
