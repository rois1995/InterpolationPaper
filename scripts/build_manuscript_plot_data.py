#!/usr/bin/env python3
"""Build manuscript figure views from archived inputs, with paired and physical-time checks."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCHEMES = ('NN1', 'BAR2', 'GP2', 'CV2')


def build(root=ROOT):
    manifest = json.loads((root/'provenance/figure_inputs.json').read_text())
    for name, digest in manifest['archived_sha256'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Figure input digest mismatch: '+name)
    def read(name):
        return pd.read_csv(root/name)
    def write(name, frame):
        if frame.empty or not np.isfinite(frame.select_dtypes(include='number')).all().all():
            raise ValueError('Empty or non-finite plot view: '+name)
        frame.to_csv(root/'generated'/name, index=False, float_format='%.12g')

    cases = read('data/riemann_lc/per_case_metrics.csv')
    selected = cases[(cases.dIters == 2) & (cases.SF == 2) & cases.scheme.isin(SCHEMES)]
    if len(selected) != 24 or selected.duplicated(['scheme','seed','lc']).any():
        raise ValueError('Expected 12 unique matched figure pairs')
    mass = []
    curves = []
    series = read('data/riemann_lc/timeseries.csv')
    for i, scheme in enumerate(SCHEMES, 1):
        paired = []
        for seed in (1, 2, 3):
            pair = selected[(selected.scheme == scheme) & (selected.seed == seed)].set_index('lc')
            standard, lc = pair.loc['off'], pair.loc['on']
            mass.append(dict(x=i+(seed-2)*0.09, scheme=scheme, seed=seed,
                transfer=abs(standard.cum_density), standard=abs(standard.mass_drift_post399),
                lc=abs(lc.mass_drift_post399)))
            def history(case):
                f = series[(series.case == case) & ((series.state == 'pre') & (series.step < 399) | (series.state == 'post') & (series.step == 399))]
                if f.duplicated('step').any():
                    raise ValueError('Ambiguous pre/post history: '+case)
                return f.set_index('step').sort_index()
            a, b = history(standard['case']), history(lc['case'])
            if not a.index.equals(b.index) or not np.array_equal(a.t, b.t) or len(a) != 42:
                raise ValueError('Unaligned paired error histories')
            change = 100*(b.L1_rho/a.L1_rho-1)
            endpoint = 100*(lc.L1_post399/standard.L1_post399-1)
            if not np.isclose(change.loc[399], endpoint, atol=0.001):
                raise ValueError('Error-curve endpoint disagrees with double-precision reduction')
            paired.append(change)
        values = pd.concat(paired, axis=1)
        curves.append(pd.DataFrame({'t':a.t, scheme+'_mean':values.mean(axis=1),
                                   scheme+'_min':values.min(axis=1), scheme+'_max':values.max(axis=1)}).reset_index(drop=True))
    mass = pd.DataFrame(mass)
    if (mass[['transfer','standard','lc']] <= 0).any().any():
        raise ValueError('Log-mass figure requires positive values; do not invent a plotting floor')
    write('riemann_mass_plot.csv',mass)
    curve = curves[0]
    for other in curves[1:]:
        if not np.array_equal(curve.t, other.t):
            raise ValueError('Unaligned scheme histories')
        curve = curve.merge(other,on='t',validate='one_to_one')
    write('riemann_error_plot.csv',curve)

    # Restart exports repeat the last two identical snapshots; reject conflicting duplicates.
    shells = read('data/euler_pressure_shells.csv').drop_duplicates()
    if shells.duplicated(['scheme','radius','step']).any():
        raise ValueError('Duplicate pressure sample')
    trace = pd.DataFrame({'lag':range(400)})
    arrival = []
    for (scheme,radius), group in shells.groupby(['scheme','radius']):
        group = group.sort_values('step')
        if not np.array_equal(group.step.to_numpy(),np.arange(400,800)) or group.pressure_rms.max() <= 0:
            raise ValueError('Incomplete or zero pressure history')
        normalized = group.pressure_rms.to_numpy()/group.pressure_rms.max()
        if scheme == 'GP2' and radius in (0.005,0.04,0.08,0.1):
            trace['r'+str(radius).replace('.','p')] = normalized
        if 0.2-radius >= 0.08:
            arrival.append(dict(scheme=scheme,travel=0.2-radius,lag=int(np.flatnonzero(normalized >= .5)[0])))
    write('euler_pressure_plot.csv',trace)
    arrival = pd.DataFrame(arrival)
    fits = read('data/naca_pressure_propagation.csv')
    for scheme, method in [('GP2','ConsGalerkinProj'),('CV2','Conservative')]:
        rows = arrival[arrival.scheme == scheme].sort_values('travel')
        slope = np.polyfit(rows.travel, rows.lag, 1)[0]
        expected = fits[(fits['case'] == f'WD_0.20-0.25_{method}_OneShot_Long_Second') & (fits.threshold == 50) & (fits.travel_exclusion == .08)]
        if len(expected) != 1 or not np.isclose(slope,expected.iloc[0].slope,rtol=.001):
            raise ValueError('Pressure onsets disagree with archived propagation fit')
        write('euler_arrival_'+scheme+'.csv', rows)

    rans = read('RANS_PostProcessing/rans_control_series.csv')
    wide = pd.DataFrame({'step':range(2,40)})
    for model in ('SA','SST'):
        for setting in (('production','halfstep','quarterstep') if model == 'SA' else ('production','halfstep')):
            for quantity in ('GP2','P'):
                f = rans[(rans.model == model) & (rans.setting == setting) & (rans.quantity == quantity) & (rans.coeff == 'CD') & rans.step.between(2,39)]
                if len(f) != 38 or f.step.duplicated().any():
                    raise ValueError('Incomplete aligned RANS window')
                if quantity == 'P':
                    subset = rans[(rans.model == model) & (rans.setting == setting) & (rans.coeff == 'CD') & rans.step.between(2,39)]
                    check = subset.pivot(index='step',columns='quantity',values='value')
                    # The archived series use seven significant figures independently per column.
                    tolerance = 1e-6*(check.NN2.abs()+check.GP2.abs()+check.P.abs())+1e-15
                    if (abs(check.P-(check.NN2-check.GP2)) > tolerance).any():
                        raise ValueError('RANS difference does not equal NN2 minus GP2')
                wide = wide.merge(f[['step','value']].rename(columns={'value':f'{model}_{quantity}_{setting}'}),on='step',validate='one_to_one')
    write('rans_drag_plot.csv',wide)
    print('Validated paired error curves, final mass diagnostics, pressure arrivals and aligned RANS histories')


if __name__ == '__main__':
    build()
