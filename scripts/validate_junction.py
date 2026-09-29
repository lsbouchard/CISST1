#!/usr/bin/env python3
"""Validate every stored junction sample and recompute one full transport case."""
import csv
import hashlib
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
from junction_transport import Junction, ThermalBath, solve, energy_grid, dilute_scattering

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/junction'


def read(name):
    with (DATA/name).open() as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key, value in row.items():
            if value == '':
                row[key] = None
            elif value in ('True', 'False'):
                row[key] = value == 'True'
            else:
                try:
                    row[key] = float(value)
                    if not np.isfinite(row[key]):
                        raise AssertionError(f'nonfinite {name}: {key}')
                except ValueError:
                    pass
    return rows


def check_observables(row):
    if row['fixed_point_residual'] >= 1e-11 or row['charge_balance_relative'] > 1e-7:
        raise AssertionError('unconverged or nonconserving transport sample')
    for key, tolerance in [('energy_balance', 1e-9), ('event_energy_balance', 1e-8), ('spin_balance', 1e-9)]:
        if abs(row[key]) > tolerance:
            raise AssertionError(f'failed {key}')
    for port in ['left', 'right']:
        if abs(row[f'extraction_Pz_{port}']) > 1+1e-9:
            raise AssertionError('unphysical extraction spin density')
        particle, spin = row[f'particle_into_{port}'], row[f'spin_sigma_into_{port}']
        ratio = row[f'net_spin_per_particle_{port}']
        if ratio is None:
            if abs(particle) >= 1e-12:
                raise AssertionError('missing finite-current spin ratio')
        elif not np.isclose(ratio*particle, spin, atol=1e-12, rtol=1e-9):
            raise AssertionError('spin/particle ratio inconsistent with fluxes')
    if row['occupied_empty_identity_error'] > 1e-10:
        raise AssertionError('spectral identity failed')
    particle_scale = max(abs(row['particle_into_left']), abs(row['particle_into_right']),
                         row['gross_extraction_left']+row['gross_injection_left'],
                         row['gross_extraction_right']+row['gross_injection_right'], 1e-15)
    if abs(row['phonon_number_collision']) > 1e-7*particle_scale:
        raise AssertionError('phonons do not conserve electron number')
    if row['density_min_eigenvalue'] < -1e-9 or row['density_max_eigenvalue'] > 1.003:
        raise AssertionError('unphysical density eigenvalue')


def main():
    summary = json.loads((DATA/'parameters.json').read_text())
    for key, name in [('solver_sha256', 'junction_transport.py'),
                      ('pauli_source_sha256', 'hahn_transport.py'),
                      ('generator_sha256', 'benchmark_junction.py')]:
        if summary['provenance'][key] != hashlib.sha256((ROOT/'scripts'/name).read_bytes()).hexdigest():
            raise AssertionError('stored data belong to different source code')
    samples = {name: read(name+'.csv') for name in ['convergence', 'controls', 'feasibility', 'variations', 'dilute_scattering']}
    counts = dict(convergence=5, controls=7, feasibility=12, variations=14, dilute_scattering=16)
    for name, count in counts.items():
        if len(samples[name]) != count:
            raise AssertionError(f'incorrect sample count: {name}')
        if name != 'dilute_scattering':
            for row in samples[name]:
                check_observables(row)
    if {(r['sites'], r['eta']) for r in samples['feasibility']} != {(s, e) for s in [5, 9, 13] for e in [0., .04, .08, .16]}:
        raise AssertionError('incomplete parameter map')
    for row in samples['feasibility']+samples['variations']:
        if not row['numerically_resolved'] or row['P_grid_difference'] >= .03 or row['current_grid_difference'] >= .002 or row['spectral_sumrule_error'] >= .003:
            raise AssertionError('unresolved published parameter point')
    if summary['resolved_map_points'] != 12 or summary['total_map_points'] != 12:
        raise AssertionError('incorrect resolution metadata')
    maximum = max(abs(r['net_spin_per_particle_right']) for r in samples['feasibility'])
    np.testing.assert_allclose(summary['max_resolved_net_Pz'], maximum, atol=0., rtol=0.)
    maximum_tested = max(abs(r['net_spin_per_particle_right']) for r in samples['feasibility']+samples['variations'])
    np.testing.assert_allclose(summary['max_tested_net_Pz'], maximum_tested, atol=0., rtol=0.)
    controls = {r['check']: r for r in samples['controls']}
    for name in ['equilibrium', 'all_SOC_zero', 'elastic_SOC_only']:
        if abs(controls[name]['spin_sigma_into_right']) > 1e-10:
            raise AssertionError('zero-signal control failed')
    base = controls['reference']
    np.testing.assert_allclose(controls['opposite_enantiomer']['spin_sigma_into_right'], -base['spin_sigma_into_right'], atol=1e-12)
    np.testing.assert_allclose(controls['particle_flow_reversed']['extraction_Pz_left'], -base['extraction_Pz_right'], atol=1e-11)
    np.testing.assert_allclose(controls['rotating_spin_basis']['spin_sigma_into_right'], base['spin_sigma_into_right'], atol=1e-12)
    np.testing.assert_allclose(summary['reference']['net_spin_per_particle_right'], samples['convergence'][2]['net_spin_per_particle_right'], atol=0., rtol=0.)
    for row in samples['dilute_scattering']:
        model = replace(Junction(), sites=5, eta=row['eta'])
        expected = dilute_scattering(model, 1.2, states=int(row['phonon_states']))
        for key, value in expected.items():
            np.testing.assert_allclose(row[key], value, atol=1e-12, rtol=1e-8)
        if row['unitarity_error'] > 1e-12:
            raise AssertionError('scattering probability lost')
    # A fresh SCBA calculation checks the serialized map independently of its cache.
    reference = next(r for r in samples['feasibility'] if r['sites'] == 5 and r['eta'] == .08)
    bath = ThermalBath.ohmic(reference['bath_spacing_used'], .5, int(reference['bath_count_used']), .15)
    recomputed = solve(replace(Junction(), sites=5), bath,
                       energy_grid(spacing=reference['energy_spacing_used']), tolerance=1e-11).observables()
    for key in ['particle_into_right', 'spin_sigma_into_right', 'extraction_Pz_right', 'net_spin_per_particle_right']:
        np.testing.assert_allclose(reference[key], recomputed[key], atol=1e-11, rtol=1e-7)
    print('Validated all finite-junction data, controls, conservation, and a fresh transport case.')


if __name__ == '__main__':
    main()
