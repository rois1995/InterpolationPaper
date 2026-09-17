#!/usr/bin/env python3
"""Guard for the two tables typed into paper.tex: every number must equal the CSV it summarizes.
    tab:disc-selected      <- data/discontinuous_boundedness.csv
    tab:riemann-stress-new <- data/riemann_frequency_full.csv, data/riemann_amplitude_global.csv
Exit status 1 on any mismatch. Standard library only."""
import csv, re, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
tex = (ROOT / 'paper.tex').read_text()
rows = lambda p: list(csv.DictReader(open(ROOT / p, newline='')))
bad = []
def block(label):
    i = tex.index(label); a = tex.index('\\begin{tabular}', i); b = tex.index('\\end{tabular}', a)
    return [l.strip() for l in tex[a:b].split('\n') if '&' in l and not l.strip().startswith(('Transfer', 'Condition'))]
def num(t):
    t = t.replace('$', '').replace('\\times10^{', 'e').replace('}', '').replace('{', '').replace(' ', '')
    return float(t)
# discontinuous table
D = {(r['method'], r['order'], r['limiter'], r['enforce_cons']): r for r in rows('data/discontinuous_boundedness.csv')}
key = {'NN2, unlimited': ('NN', 'Second', 'none', 'False'), 'NN2, BJ': ('NN', 'Second', 'BJ', 'False'), 'NN2, BJ+GC': ('NN', 'Second', 'BJ', 'Global'),
       'BAR3, unlimited': ('Barycentric', 'Third', 'none', 'False'), 'BAR3, BJ+GC': ('Barycentric', 'Third', 'BJ', 'Global'),
       'GP2, hard limited': ('ConsGalerkinProj', 'Second', 'hard', 'False'), 'GP3, hard limited, vertex output': ('ConsGalerkinProj', 'Third', 'hard', 'False')}
for line in block('label{tab:disc-selected}'):
    cells = [c.strip().rstrip('\\').strip() for c in line.split('&')]
    if cells[0] not in key: bad.append(f'disc-selected: unknown row {cells[0]}'); continue
    r = D[key[cells[0]]]
    for got, want, name in ((num(cells[1]), float(r['overshoot']), 'overshoot'), (num(cells[2]), float(r['undershoot']), 'undershoot'), (num(cells[3]), float(r['cons_defect']), 'defect')):
        if abs(got - want) > 0.006 * max(abs(want), 1e-14):
            bad.append(f'disc-selected {cells[0]} {name}: typed {got} vs csv {want}')
# riemann stress table
F = {(r['scheme'], r['dIters']): r for r in rows('data/riemann_frequency_full.csv')}
A = {(r['scheme'], r['global_correction'], r['safe_factor']): r for r in rows('data/riemann_amplitude_global.csv')}
def val(scheme, cond):
    if cond == '80': return F[(scheme, '10')]
    if cond == 'sf4': return A[(scheme, 'off' if scheme != 'GP2' else 'intrinsic', '4')]
    return A[(scheme, 'off' if scheme != 'GP2' else 'intrinsic', '2')]
for line in block('label{tab:riemann-stress-new}'):
    cells = [c.strip().rstrip('\\').strip() for c in line.split('&')]
    cond = '80' if cells[0].startswith('80') else ('sf4' if 'SF4' in cells[0] else 'sf2')
    checks = [('NN1', 1), ('BAR2', 3), ('GP2', 4)]
    for scheme, k in checks:
        want = float(val(scheme, cond)['width_change_cells']); got = num(cells[k])
        if abs(got - want) > 0.006 * max(abs(want), 0.01): bad.append(f'riemann-stress {cells[0]} {scheme} width: typed {got} vs csv {want}')
    if cond != '80':
        want = float(A[('NN1', 'on', '4' if cond == 'sf4' else '2')]['width_change_cells']); got = num(cells[2])
        if abs(got - want) > 0.006 * abs(want): bad.append(f'riemann-stress {cells[0]} NN1+GC width: typed {got} vs csv {want}')
    pos = [abs(float(val(s, cond)['position_error_cells'])) for s in ('NN1', 'BAR2', 'GP2')] + ([abs(float(A[('NN1', 'on', '4' if cond == 'sf4' else '2')]['position_error_cells']))] if cond != '80' else [])
    got = num(cells[5]); want = max(pos)
    if abs(got - want) > 0.006 * max(want, 0.01): bad.append(f'riemann-stress {cells[0]} max position: typed {got} vs csv {want}')
print('\n'.join(bad) if bad else 'typed tables agree with the CSV files'); sys.exit(1 if bad else 0)
