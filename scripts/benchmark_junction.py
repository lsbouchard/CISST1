#!/usr/bin/env python3
"""Regenerate the finite-junction benchmark, controls and conditional feasibility map."""
from dataclasses import asdict, replace
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from junction_transport import (Junction, ThermalBath, solve, energy_grid,
                               dilute_scattering, adjoint)
from reciprocal_device import ROOT, METADATA, plt

DATA = ROOT/'data/junction'
FIGURES = ROOT/'figures/junction'


def write_rows(name, rows):
    with (DATA/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def elastic_dwell(model, energy=1.):
    """Spin-averaged elastic injectivity time [hbar/t], not a directed active time."""
    H, _ = model.operators()
    sigma, gamma = model.embeddings(np.array([energy]))
    GR = np.linalg.inv(energy*np.eye(len(H))-H-sigma[:, 0].sum(axis=0))
    return float(np.trace(GR@gamma[0, 0]@adjoint(GR)).real/2)


def main(resume=False):
    DATA.mkdir(parents=True, exist_ok=True)
    # Partial reruns must not retain a completed metadata record from older data.
    (DATA/'parameters.json').unlink(missing_ok=True)
    base = Junction()
    cache = {}
    checkpoint = ROOT/'build/junction_checkpoint.json'
    fingerprint = dict(solver_sha256=hashlib.sha256((ROOT/'scripts/junction_transport.py').read_bytes()).hexdigest(),
                       pauli_source_sha256=hashlib.sha256((ROOT/'scripts/hahn_transport.py').read_bytes()).hexdigest(),
                       numpy=np.__version__, scipy=scipy.__version__, python=platform.python_version(),
                       tolerance=1e-11, mixing=.35, include_hartree=True)
    if resume and checkpoint.exists():
        previous = json.loads(checkpoint.read_text())
        if previous['fingerprint'] == fingerprint:
            cache = previous['results']
        else:
            print('Ignoring checkpoint from a different solver/runtime.', flush=True)
    def calculate(model=base, temperature=.15, potentials=(1.4, .6), spacing=.01,
                  bath_spacing=.05, bath_count=60, window=(-5., 9.), basis=None):
        key = json.dumps([asdict(model), temperature, potentials, spacing, bath_spacing, bath_count, window], sort_keys=True)
        if basis is None and key in cache:
            return cache[key]
        start = time.perf_counter()
        bath = ThermalBath.ohmic(bath_spacing, .5, bath_count, temperature)
        result = solve(model, bath, energy_grid(*window, spacing), potentials,
                       tolerance=1e-11, basis=basis)
        obs = result.observables()
        if obs['charge_balance_relative'] > 1e-7 or abs(obs['event_energy_balance']) > 1e-8 or abs(obs['spin_balance']) > 1e-9:
            raise AssertionError('conservation check failed')
        print(f'sites={model.sites}, eta={model.eta:g}, zeta={model.zeta:g}, contacts={model.contact_left:g}/{model.contact_right:g}, dE={spacing:g}: Pnet,R={obs["net_spin_per_particle_right"]}, {time.perf_counter()-start:.1f}s', flush=True)
        if basis is None:
            cache[key] = obs
            checkpoint.parent.mkdir(exist_ok=True)
            temporary = checkpoint.with_suffix('.tmp')
            temporary.write_text(json.dumps(dict(fingerprint=fingerprint, results=cache), indent=2, allow_nan=False)+'\n')
            temporary.replace(checkpoint)
        return obs

    convergence = []
    for label, spacing, bath_spacing, count, window in [
        ('coarse', .02, .1, 30, (-5., 9.)),
        ('medium', .01, .05, 60, (-5., 9.)),
        ('fine', .005, .025, 120, (-5., 9.)),
        ('bath_tail', .005, .025, 160, (-5., 9.)),
        ('energy_window', .005, .025, 120, (-7., 11.))]:
        obs = calculate(spacing=spacing, bath_spacing=bath_spacing, bath_count=count, window=window)
        convergence.append(dict(check=label, energy_spacing=spacing, bath_spacing=bath_spacing,
                                bath_max_energy=count*bath_spacing, energy_low=window[0], energy_high=window[1], **obs))
        write_rows('convergence.csv', convergence)
    reference = convergence[2]
    fine_P = reference['net_spin_per_particle_right']
    differences = {r['check']: abs(r['net_spin_per_particle_right']-fine_P)/abs(fine_P)
                   for r in convergence if r['check'] != 'fine'}
    if max(differences[name] for name in ['medium', 'bath_tail', 'energy_window']) > .03:
        raise AssertionError(f'reference discretization not converged to 3%: {differences}')

    controls = []
    for label, model, potentials, basis in [
        ('reference', base, (1.4, .6), None),
        ('equilibrium', base, (1., 1.), None),
        ('all_SOC_zero', replace(base, zeta=0., eta=0.), (1.4, .6), None),
        ('elastic_SOC_only', replace(base, eta=0.), (1.4, .6), None),
        ('opposite_enantiomer', replace(base, chi=-1), (1.4, .6), None),
        ('particle_flow_reversed', base, (.6, 1.4), None)]:
        controls.append(dict(check=label, energy_spacing=.01, bath_spacing=.05,
                             bath_count=60, temperature=.15, mu_left=potentials[0],
                             mu_right=potentials[1], **calculate(model, potentials=potentials)))
    phi = base.chi*(np.arange(base.sites)-(base.sites-1)/2)*base.step
    U = np.diag(np.column_stack([np.exp(-.5j*phi), np.exp(.5j*phi)]).ravel())
    controls.append(dict(check='rotating_spin_basis', energy_spacing=.01,
                         bath_spacing=.05, bath_count=60, temperature=.15,
                         mu_left=1.4, mu_right=.6, **calculate(basis=U)))
    write_rows('controls.csv', controls)

    def checked(model=base, temperature=.15, potentials=(1.4, .6)):
        coarse = calculate(model, temperature, potentials, spacing=.02, bath_spacing=.1, bath_count=30)
        medium = calculate(model, temperature, potentials)
        selected_spacing, selected_bath_spacing, selected_bath_count = .01, .05, 60
        current_error = abs(medium['particle_into_right']-coarse['particle_into_right'])/max(abs(medium['particle_into_right']), 1e-12)
        P = medium['net_spin_per_particle_right']
        polarization_error = abs(P-coarse['net_spin_per_particle_right'])/max(abs(P), 1e-10)
        if polarization_error > .03 or current_error > .002 or medium['spectral_sumrule_error'] > .003:
            previous = medium
            medium = calculate(model, temperature, potentials, spacing=.005, bath_spacing=.025, bath_count=120)
            selected_spacing, selected_bath_spacing, selected_bath_count = .005, .025, 120
            current_error = abs(medium['particle_into_right']-previous['particle_into_right'])/max(abs(medium['particle_into_right']), 1e-12)
            P = medium['net_spin_per_particle_right']
            polarization_error = abs(P-previous['net_spin_per_particle_right'])/max(abs(P), 1e-10)
        valid = polarization_error < .03 and current_error < .002 and medium['spectral_sumrule_error'] < .003
        if not valid:
            raise AssertionError(f'unresolved map/variation: {model}, T={temperature}, mu={potentials}')
        return dict(current_grid_difference=current_error, P_grid_difference=polarization_error,
                    numerically_resolved=valid, energy_spacing_used=selected_spacing,
                    bath_spacing_used=selected_bath_spacing, bath_count_used=selected_bath_count, **medium)

    map_rows = []
    for sites in [5, 9, 13]:
        for eta in [0., .04, .08, .16]:
            model = replace(base, sites=sites, eta=eta)
            map_rows.append(dict(sites=sites, geometrical_turns=model.turns, eta=eta,
                                 elastic_reference_dwell=elastic_dwell(model),
                                 **checked(model)))
            write_rows('feasibility.csv', map_rows)

    variations = []
    for name, value, model, temperature, potentials in [
        ('contact', .7, replace(base, contact_left=.7, contact_right=.7), .15, (1.4, .6)),
        ('contact', 1., base, .15, (1.4, .6)),
        ('contact', 1.4, replace(base, contact_left=1.4, contact_right=1.4), .15, (1.4, .6)),
        ('temperature', .075, base, .075, (1.4, .6)),
        ('temperature', .15, base, .15, (1.4, .6)),
        ('temperature', .3, base, .3, (1.4, .6)),
        ('SOC_scale', 0., replace(base, zeta=0., eta=0.), .15, (1.4, .6)),
        ('SOC_scale', .5, replace(base, zeta=.075, eta=.04), .15, (1.4, .6)),
        ('SOC_scale', 1., base, .15, (1.4, .6)),
        ('SOC_scale', 2., replace(base, zeta=.3, eta=.16), .15, (1.4, .6)),
        ('bias', .2, base, .15, (1.1, .9)),
        ('bias', .8, base, .15, (1.4, .6)),
        ('bias', 1.6, base, .15, (1.8, .2)),
        ('dilute_drain', -3., base, .15, (1.4, -3.))]:
        obs = checked(model, temperature, potentials)
        variations.append(dict(parameter=name, value=value, temperature=temperature,
                               mu_left=potentials[0], mu_right=potentials[1],
                               elastic_reference_dwell=elastic_dwell(model), **obs))
        write_rows('variations.csv', variations)

    scattering = []
    for states in [2, 3, 4, 6]:
        for eta in [0., .01, .02, .04]:
            model = replace(base, sites=5, eta=eta)
            scattering.append(dict(phonon_states=states, eta=eta, **dilute_scattering(model, 1.2, states=states)))
    write_rows('dilute_scattering.csv', scattering)
    resolved = [r for r in map_rows if r['numerically_resolved']]
    summary = dict(status='synthetic finite-junction benchmark; not an experimental fit or material identification',
                   method='Full spin-matrix Fock+Hartree SCBA; spin-degenerate finite-band Fermi leads; thermal Ohmic amplitude bath',
                   provenance={**fingerprint, 'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                   units=dict(energy='t', time='hbar/t', particle_current='t/hbar', spin_rate='sigma_z-weighted particle current; multiply by hbar/2',
                              spin_torque='t', heat_current='t^2/hbar'),
                   reference_model=asdict(base), reference_temperature=.15, reference_mu=[1.4, .6],
                   elastic_reference_dwell_energy=1.,
                   reference=reference, reference_relative_discretization_differences=differences,
                   max_resolved_net_Pz=float(max(abs(r['net_spin_per_particle_right']) for r in resolved)),
                   max_tested_net_Pz=float(max(abs(r['net_spin_per_particle_right']) for r in map_rows+variations)),
                   resolved_map_points=len(resolved), total_map_points=len(map_rows),
                   scientific_limits=['SCBA is a weak-coupling approximation; numerical convergence does not certify its physical accuracy',
                                      'No unique scalar T1 is extracted from a steady-state calculation',
                                      'Gross extraction polarization is not the same as source-resolved transmission',
                                      'Elastic reference dwell time is not a directed active alignment time',
                                      'Material SOC, effective mass and phonon vertex remain unmeasured inputs'])
    (DATA/'parameters.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    write_figures(map_rows, variations, controls)
    print(json.dumps(summary, indent=2, allow_nan=False), flush=True)


def write_figures(map_rows, variations, controls):
    FIGURES.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.9), layout='constrained')
    for eta, color in [(.04, '#287c83'), (.08, '#ad3d38'), (.16, '#58559a')]:
        rows = [r for r in map_rows if r['eta'] == eta]
        axes[0].plot([r['geometrical_turns'] for r in rows], [100*r['net_spin_per_particle_right'] for r in rows],
                     marker='o', color=color, label=f'$\\eta={eta:g}$')
    axes[0].set(xlabel='Geometrical helical turns', ylabel='Net drain spin / particle (%)', title='(a) Finite-junction response')
    axes[0].legend(frameon=False, fontsize=8)
    contact = [r for r in variations if r['parameter'] == 'contact']
    axes[1].plot([r['value'] for r in contact], [100*r['net_spin_per_particle_right'] for r in contact], color='#287c83', marker='o', label='Net drain ratio')
    axes[1].plot([r['value'] for r in contact], [100*r['extraction_Pz_right'] for r in contact], color='#ad3d38', marker='s', label='Gross extraction')
    axes[1].set(xlabel='Contact hopping / t', ylabel='Physical-spin signal (%)', title='(b) Contacts and observable')
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES/'fig_junction_feasibility.pdf', metadata={**METADATA, 'Title': 'Conditional finite-junction feasibility'})
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.9), layout='constrained')
    controls_by_name = {r['check']: r for r in controls}
    checks = [controls_by_name[name] for name in
              ['reference', 'opposite_enantiomer', 'all_SOC_zero', 'elastic_SOC_only']]
    labels = ['Reference', r'$\chi=-1$', 'Zero SOC', 'Elastic']
    axes[0].bar(np.arange(len(checks)), [100*r['net_spin_per_particle_right'] for r in checks], color=['#287c83', '#ad3d38', '#737373', '#737373'])
    axes[0].axhline(0, color='.4', lw=.6)
    axes[0].set(xticks=np.arange(len(checks)), xticklabels=labels, ylabel='Net drain spin / particle (%)', title='(a) Symmetry and bath controls')
    axes[0].tick_params(axis='x', labelsize=8)
    for parameter, color, symbol in [('temperature', '#287c83', 'o'), ('SOC_scale', '#ad3d38', 's')]:
        rows = [r for r in variations if r['parameter'] == parameter]
        x = [r['value']/.15 if parameter == 'temperature' else r['value'] for r in rows]
        axes[1].plot(x, [100*r['net_spin_per_particle_right'] for r in rows], color=color, marker=symbol,
                     label='Temperature / reference' if parameter == 'temperature' else 'Joint SOC scale')
    axes[1].set(xlabel='Parameter / reference value', ylabel='Net drain spin / particle (%)', title='(b) SOC and temperature')
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES/'fig_junction_controls.pdf', metadata={**METADATA, 'Title': 'Finite-junction controls and parameter dependence'})
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--resume', action='store_true', help='reuse checkpoints only when solver, runtime and inputs match')
    modes.add_argument('--plots-only', action='store_true', help='validate stored data and regenerate figures without repeating the full scan')
    args = parser.parse_args()
    if args.plots_only:
        from validate_junction import main as validate, read
        validate()
        write_figures(read('feasibility.csv'), read('variations.csv'), read('controls.csv'))
    else:
        main(resume=args.resume)
