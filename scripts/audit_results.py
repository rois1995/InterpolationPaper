#!/usr/bin/env python3
"""Audit the tabulated evidence without modifying any input or generated CSV.

This script is not a CFD validation run. It recomputes quantities from the
available summary values, flags invalid reported diagnostics, and checks the
manuscript's citation/label consistency. Python 3.10+; standard library only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys
from typing import Any


def load_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required audit input is missing: {path}")
    with path.open(newline='', encoding='utf-8-sig') as stream:
        result = list(csv.DictReader(stream))
    if not result:
        raise ValueError(f"No data rows in {path}")
    return result


def linear_fit(x: list[float], y: list[float]) -> dict[str, float]:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("A fit needs at least two paired samples")
    if not all(math.isfinite(v) for v in x + y):
        raise ValueError("Non-finite fit input")
    xb, yb = statistics.mean(x), statistics.mean(y)
    denominator = sum((v - xb)**2 for v in x)
    if denominator == 0:
        raise ValueError("Cannot fit coincident abscissae")
    slope = sum((a-xb)*(b-yb) for a, b in zip(x, y)) / denominator
    intercept = yb - slope*xb
    ss = sum((b - slope*a - intercept)**2 for a, b in zip(x, y))
    total = sum((b-yb)**2 for b in y)
    return {"slope": slope, "intercept": intercept,
            "R_squared": 1 - ss/total if total else 1.0}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, default=Path('audit/numerical_audit.json'))
    parser.add_argument('--check-preservation', action='store_true',
                        help='Fail if a reviewed CSV differs from the original snapshot')
    args = parser.parse_args()
    root = args.root.resolve()
    read = lambda name: load_rows(root / name)
    out: dict[str, Any] = {}
    out['scope'] = 'Recomputation from archived tables; no raw-field/solver validation.'
    out['csv_sha256'] = {str(f.relative_to(root)): hashlib.sha256(f.read_bytes()).hexdigest()
                         for folder in ('data', 'generated')
                         for f in sorted((root/folder).glob('*.csv'))}
    manifest = root/'audit/original_csv_sha256.json'
    if manifest.exists():
        original = json.loads(manifest.read_text())
        out['changed_or_missing_original_csvs'] = [k for k,v in original.items()
                                                   if out['csv_sha256'].get(k) != v]
    else:
        out['changed_or_missing_original_csvs'] = ['manifest_missing']

    rows = read('data/naca_global_correction.csv')
    out['naca_gc_relative_changes_percent'] = {}
    for raw, repaired in zip(rows[::2], rows[1::2]):
        out['naca_gc_relative_changes_percent'][raw['method']] = {
            name: 100*(float(repaired[name])-float(raw[name]))/float(raw[name])
            for name in ('epsilon_norm','delta_full_norm','delta_CL_full')}

    rows = read('data/naca_acoustic_delay.csv')
    x = [float(row['distance']) for row in rows]
    out['force_onset_fits'] = {k: linear_fit(x,[float(row[k]) for row in rows])
                               for k in ('onset_1','onset_5','onset_25','onset_50')}
    rows = read('data/naca_pressure_propagation.csv')
    slopes = [float(row['slope']) for row in rows]
    speeds = [1/(0.001*s) for s in slopes]
    out['pressure_threshold_sensitivity'] = {
        'mean_slope':statistics.mean(slopes), 'slope_sample_sd':statistics.stdev(slopes),
        'mean_speed':statistics.mean(speeds), 'speed_sample_sd':statistics.stdev(speeds),
        'speed_min':min(speeds), 'speed_max':max(speeds),
        'interpretation':'Six thresholds of one event, not six independent repetitions.'}
    rows = read('data/naca_force_drift.csv')
    means = {m:statistics.mean(float(r['CD_trend']) for r in rows if r['method']==m)
             for m in ('NN2','GP2')}
    out['drag_trend_means'] = means
    out['drag_trend_ratio_NN2_GP2'] = means['NN2']/means['GP2']
    rows = read('data/naca_initial_residual_identity.csv')
    out['invalid_cosine_rows'] = [r for r in rows
                                  if not (-1 <= float(r['cosine_min']) <= 1)
                                  or not (-1 <= float(r['cosine_max']) <= 1)]
    rows = read('data/naca_transfer_defects.csv')
    ranked = sorted((r['method'],float(r['epsilon_norm'])) for r in rows)
    out['naca_epsilon_values'] = dict(ranked)
    values = [v for m,v in ranked if m!='NN1']
    out['naca_higher_order_spread_percent'] = 100*(max(values)/min(values)-1)
    rows = read('data/naca_force_history.csv')
    out['relaxation_summaries'] = {
        r['method']:{'peak_over_first':float(r['max_abs_delta_CL'])/abs(float(r['delta_CL_first'])),
                     'last_abs_over_first_abs':abs(float(r['delta_CL_last'])/float(r['delta_CL_first']))}
        for r in rows}
    rows = read('data/naca_first_step_response.csv')
    out['first_step_ratios'] = {}
    for method in dict.fromkeys(r['method'] for r in rows):
        bycase = {r['case']:r for r in rows if r['method']==method}
        out['first_step_ratios'][method] = {
            name:{case:float(bycase[case][name])/float(bycase['full'][name])
                  for case in ('latest','older')}
            for name in ('delta_norm','delta_CL')}
    rows = read('data/riemann_refinement_errors.csv')
    x = [math.log(float(r['h'])) for r in rows]
    out['riemann_order_fits'] = {m:linear_fit(x,[math.log(float(r[m])) for r in rows])
                                for m in rows[0] if m!='h'}
    rows = read('data/riemann_global_correction.csv')
    out['riemann_gc_error_changes_percent'] = {
        r['scheme']:100*(float(r['l1_global'])/float(r['l1_off'])-1) for r in rows}
    out['static_orders'] = {name:read(f'generated/{name}_orders.csv')
                           for name in ('asymmetric','discontinuous')}
    # Independent finite-sum check of the manuscript's mathematical counterexample.
    out['quadratic_lumping_counterexample'] = []
    for n in (4,8,16,32):
        h=1/n
        q=h*(0.5+sum((i*h)**2 for i in range(1,n)))
        out['quadratic_lumping_counterexample'].append(
            {'n':n,'computed_sum':q,'error':q-1/3,'formula_error':h*h/6})
        if not math.isclose(q-1/3,h*h/6,rel_tol=1e-11,abs_tol=1e-15):
            raise AssertionError('Analytical counterexample check failed')

    tex=(root/'paper.tex').read_text()
    bib=(root/'References.bib').read_text()
    keys=re.findall(r'@\w+\s*\{\s*([^,]+),',bib)
    used=[]
    for group in re.findall(r'\\(?:cite|parencite|textcite|autocite)(?:\[[^\]]*\])*\{([^}]+)\}',tex):
        used += [s.strip() for s in group.split(',')]
    labels=re.findall(r'\\label\{([^}]+)\}',tex)
    refs=re.findall(r'\\(?:ref|eqref)\{([^}]+)\}',tex)
    # pgfplots named legends are valid LaTeX references without an explicit \label.
    legends=re.findall(r'legend to name=([^,\]\s]+)',tex)
    out['source_checks']={
        'bibliography_entries':len(keys),'cited_keys':len(set(used)),
        'duplicate_bibliography_keys':sorted({k for k in keys if keys.count(k)>1}),
        'missing_citation_keys':sorted(set(used)-set(keys)),
        'duplicate_labels':sorted({k for k in labels if labels.count(k)>1}),
        'undefined_reference_labels':sorted(set(refs)-set(labels)-set(legends))}
    output=args.output if args.output.is_absolute() else root/args.output
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(f'Audited {len(out["csv_sha256"])} CSVs; report: {output}')
    print(f'Invalid cosine summary rows retained in source: {len(out["invalid_cosine_rows"])}')
    print('Source checks:',out['source_checks'])
    failed=any(out['source_checks'][name] for name in (
        'duplicate_bibliography_keys','missing_citation_keys','duplicate_labels',
        'undefined_reference_labels'))
    if args.check_preservation and out['changed_or_missing_original_csvs']:
        print('Preservation mismatch:',out['changed_or_missing_original_csvs'],file=sys.stderr)
        failed=True
    return 1 if failed else 0


if __name__=='__main__':
    try:
        raise SystemExit(main())
    except (OSError,ValueError,KeyError,AssertionError) as exc:
        print(f'Audit failed: {exc}',file=sys.stderr)
        raise SystemExit(2)
