#!/usr/bin/env python3
"""Build Riemann flux-correction tables from archived paired results; never edit inputs."""
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/riemann_lc'

def rows(name):
    with (SOURCE/name).open(newline='') as stream:
        return list(csv.DictReader(stream))

def write(name, records):
    with (ROOT/'generated'/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

def main():
    manifest = json.loads((ROOT/'provenance/riemann_lc_inputs.json').read_text())
    for name, digest in manifest['sha256'].items():
        if hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Riemann input digest mismatch: '+name)
    cases = rows('per_case_metrics.csv')
    pairs = rows('paired_lc_vs_off.csv')
    keyed = {(r['twin'], r['lc']): r for r in cases}
    if len(keyed) != 72 or len(cases) != 72 or len(pairs) != 36 or len({r['twin'] for r in pairs}) != 36:
        raise ValueError('Expected 36 unique matched Riemann pairs')
    for pair in pairs:
        standard, lc = (keyed[(pair['twin'], flag)] for flag in ('off','on'))
        for key in ('scheme','dIters','SF','seed'):
            if standard[key] != lc[key] or standard[key] != pair[key]:
                raise ValueError('Unmatched Riemann configuration: '+pair['twin'])
        for row in (standard,lc):
            if not all(math.isfinite(float(row[k])) for k in ('wid399_narrow','mass_drift_post399','energy_drift_post399','L1_post399')):
                raise ValueError('Non-finite Riemann result')
        change = 100*(float(lc['L1_post399'])/float(standard['L1_post399'])-1)
        if not math.isclose(change,float(pair['L1_post399_change_pct']),abs_tol=0.001):
            raise ValueError('Inconsistent paired L1 change: '+pair['twin'])
    summary=[]
    for scheme in ('NN1','BAR2','GP2','CV2'):
        for sf in ('4','2'):
            selected=[r for r in pairs if r['scheme']==scheme and r['SF']==sf and r['dIters']=='2']
            if {r['seed'] for r in selected} != {'1','2','3'} or len(selected)!=3:
                raise ValueError('Expected three seeds for '+scheme+' s='+sf)
            summary.append(dict(scheme=scheme,s=sf,
                width_standard=statistics.mean(float(r['wid399_narrow_off']) for r in selected),
                width_lc=statistics.mean(float(r['wid399_narrow_LC']) for r in selected),
                mass_standard=max(abs(float(r['mass_drift_post399_off'])) for r in selected),
                mass_lc=max(abs(float(r['mass_drift_post399_LC'])) for r in selected),
                L1_change=statistics.mean(float(r['L1_post399_change_pct']) for r in selected)))
    write('riemann_lc_summary.csv',summary)
    detailed=[]
    for r in sorted(pairs,key=lambda r:(-int(r['dIters']),-int(r['SF']),r['scheme'],int(r['seed']))):
        detailed.append(dict(scheme=r['scheme'],interval=r['dIters'],s=r['SF'],seed=r['seed'],
            L1_change=r['L1_post399_change_pct'],position_standard=r['pos399_off'],position_lc=r['pos399_LC'],
            width_standard=r['wid399_narrow_off'],width_lc=r['wid399_narrow_LC']))
    write('riemann_lc_pairs.csv',detailed)
    budgets={r['case']:r for r in rows('mass_budget_step399.csv')}
    sample=[]
    for scheme,sf in [('NN1','4'),('GP2','2'),('CV2','2')]:
        for flag,label in [('off','Standard'),('on','LC')]:
            chosen=[r for r in cases if r['scheme']==scheme and r['SF']==sf and r['seed']=='1' and r['dIters']=='2' and r['lc']==flag]
            if len(chosen)!=1:raise ValueError('Ambiguous mass budget case')
            r=budgets[chosen[0]['case']]
            sample.append(dict(scheme=scheme,s=sf,treatment=label,cfd_step=r['cfd_step_399'],transfer=r['last_transfer'],drift=r['total_drift']))
    write('riemann_lc_budget.csv',sample)
    print('Validated 36 matched pairs; built shock-only, flux-correction and mass-budget tables')

if __name__ == '__main__':
    main()
