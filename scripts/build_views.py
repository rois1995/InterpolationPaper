#!/usr/bin/env python3
"""Manuscript views of the RANS records (generated/reassessed_*.csv), derived from RANS_PostProcessing/*.csv.

    reassessed_rans_integrals.csv  maximum |source-to-target integral defect| at a transfer event of the
                                   repeated-remeshing runs, per closure and inventory, NN2 against GP2
    reassessed_rans_inner.csv      inner-iteration record of the repeated-remeshing runs
    reassessed_rans_controls.csv   change of the case-minus-reference response and of the paired
                                   NN2-GP2 difference between the production setting and each control
    reassessed_rans_finemesh.csv   no-transfer steady loads on the production and the finer mesh
Standard library only. Run from the repository root (the Makefile does).
"""
import csv
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'RANS_PostProcessing')
GEN = os.path.join(ROOT, 'generated')


def rows(name):
    p = os.path.join(SRC, name)
    return list(csv.DictReader(open(p, newline=''))) if os.path.exists(p) else []


def write(name, header, data):
    os.makedirs(GEN, exist_ok=True)
    p = os.path.join(GEN, name)
    text = ','.join(header) + '\n' + ''.join(','.join(str(v) for v in r) + '\n' for r in data)
    if not os.path.exists(p) or open(p).read() != text:
        open(p, 'w').write(text)
    print(f'{name}: {len(data)} rows')


def integrals():
    mx = defaultdict(float)
    for x in rows('rans_transfer_integrals.csv'):
        if 'Drift' not in x['case']:
            continue
        m = 'GP2' if 'ConsGalerkinProj' in x['case'] else 'NN2'
        try:
            v = abs(float(x['rel_defect_(tgt-src)/src']))
        except ValueError:
            continue
        mx[(x['model'], x['field'], m)] = max(mx[(x['model'], x['field'], m)], v)
    label = {'Density': 'rho', 'Momentum_x': 'rho u', 'Momentum_y': 'rho v', 'Energy': 'rho E', 'Nu_Tilde': 'rho nu-tilde', 'Turb_Kin_Energy': 'rho k', 'Omega': 'rho omega'}
    out = []
    for model in ('SA', 'SST'):
        for f in ('Density', 'Momentum_x', 'Momentum_y', 'Energy', 'Nu_Tilde', 'Turb_Kin_Energy', 'Omega'):
            a, b = mx.get((model, f, 'NN2')), mx.get((model, f, 'GP2'))
            if a is not None and b:
                out.append([model, label[f], f'{a:.3e}', f'{b:.3e}', f'{a / b:.1f}'])
    write('reassessed_rans_integrals.csv', ['closure', 'quantity', 'NN2_max', 'GP2_max', 'ratio'], out)


def inner():
    out = []
    for x in rows('rans_inner_convergence.csv'):
        if 'Drift' not in x['run']:
            continue
        out.append([x['run'].split('_')[1], 'GP2' if 'ConsGalerkinProj' in x['run'] else 'NN2', 2 if 'Seed2' in x['run'] else 1,
                    x['physical_steps'], x['steps_at_inner_cap'], x['end_median'], x['end_max'], x['INNER_ITER'], x['CONV_RESIDUAL_MINVAL']])
    write('reassessed_rans_inner.csv', ['closure', 'method', 'seed', 'steps', 'steps_at_cap', 'end_residual_median', 'end_residual_worst', 'inner_cap', 'target_log10_rms_rho'], out)


def controls():
    out = [[x['model'], x['setting'].replace(' - production', ''), x['quantity'], x['coeff'], x['rms_of_production_signal'], x['rms'], x['rms_change_over_production_rms']]
           for x in rows('rans_control_sensitivity.csv') if x['setting'].endswith('- production') and x['coeff'] in ('CL', 'CD', 'CDv')]
    if out:
        write('reassessed_rans_controls.csv', ['closure', 'control', 'quantity', 'coeff', 'production_rms', 'change_rms', 'change_over_production'], out)


def finemesh():
    out = [[x['model'], x['coeff'], x['coarse'], x['fine'], x['fine_minus_coarse'], x['percent']] for x in rows('rans_fine_mesh.csv') if x['quantity'] == 'steady']
    if out:
        write('reassessed_rans_finemesh.csv', ['closure', 'coeff', 'coarse', 'fine', 'difference', 'percent'], out)


if __name__ == '__main__':
    integrals(); inner(); controls(); finemesh()
