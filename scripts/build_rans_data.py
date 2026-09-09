#!/usr/bin/env python3
"""Audit RANS exports and prepare CSV-driven LaTeX views without modifying inputs.

The surface exports repeat identical snapshots. One record per (case, physical step)
is used for statistics; conflicting duplicates are an error. Original CSVs are retained.
Unadjusted OLS errors are descriptive fit diagnostics, not ensemble uncertainties.
The integral diagnostic is read as exported: its small values cannot be reconstructed
by subtracting the separately rounded, large source/target integrals.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import re
from pathlib import Path
import numpy as np
import pandas as pd


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text() != text:
        path.write_text(text)


def write_csv(path: Path, frame: pd.DataFrame) -> None:
    write_text(path, frame.to_csv(index=False, float_format='%.12e', na_rep='nan', lineterminator='\n'))


def case_info(case: str) -> tuple[str, str, int]:
    m = re.match(r'RANS_(SA|SST)_WD_', case)
    if m is None:
        raise ValueError(f'Unrecognized case: {case}')
    method = 'GP2' if '_ConsGalerkinProj_' in case else 'NN2' if '_NN_' in case else None
    if method is None or not case.endswith('_Second'):
        raise ValueError(f'Unexpected method/order: {case}')
    return m.group(1), method, 2 if '_Seed2_' in case else 1


def ols(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 3 or np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        raise ValueError('OLS requires at least three finite observations')
    xc = x-x.mean()
    slope = float(xc @ (y-y.mean())/(xc @ xc))
    residual = y-(y.mean()+slope*xc)
    se = float(np.sqrt((residual @ residual)/(len(x)-2)/(xc @ xc)))
    corr = float(np.corrcoef(residual[:-1], residual[1:])[0,1]) if np.std(residual)>0 else 0.0
    return dict(trend=100*slope, se=100*se, mean=float(y.mean()), rms=float(np.sqrt(np.mean(y*y))),
                peak=float(np.max(np.abs(y))), lag1_corr=corr)


def first_number(value: object) -> float:
    match = re.match(r'^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)',str(value))
    if not match:
        raise ValueError(f'Expected a number at start of {value!r}')
    return float(match.group(1))


def read_markdown_table(text: str, header: str) -> pd.DataFrame:
    lines = text.splitlines()
    start = next(i for i,l in enumerate(lines) if l.startswith(header))
    cols = [x.strip() for x in lines[start].strip().strip('|').split('|')]
    rows=[]
    for l in lines[start+2:]:
        if not l.startswith('|'): break
        row=[x.strip() for x in l.strip().strip('|').split('|')]
        if len(row)!=len(cols): raise ValueError('Malformed report table')
        rows.append(row)
    return pd.DataFrame(rows,columns=cols)


def build(root: Path) -> None:
    src=root/'RANS_PostProcessing'; gen=root/'generated'; audit_dir=root/'audit'
    frames:dict[str,pd.DataFrame]={}; audits={}; hashes={}
    for path in sorted(src.glob('*.csv')):
        frame=pd.read_csv(path)
        n_raw=len(frame)
        frame=frame.drop_duplicates().copy()
        key=['case','step'] if {'case','step'}.issubset(frame.columns) else None
        if path.name=='rans_surface_profiles.csv': key=['case','step','x','y']
        if path.name=='rans_transfer_integrals.csv': key=['case','step','which','field']
        if key is not None and frame.duplicated(key).any():
            raise ValueError(f'Conflicting duplicate physical keys in {path.name}')
        frames[path.name]=frame
        audits[path.name]={'input_rows':n_raw,'unique_rows':len(frame),'exact_duplicates_removed_in_view':n_raw-len(frame)}
        hashes[str(path.relative_to(root))]=hashlib.sha256(path.read_bytes()).hexdigest()
    required=['rans_drift_series.csv','rans_surface_forces.csv','rans_transfer_integrals.csv','rans_oneshot_lags.csv',
              'rans_pressure_probe.csv','rans_turb_state_norms_SA.csv','rans_turb_state_norms_SST.csv']
    if not all(n in frames for n in required): raise FileNotFoundError('Missing required RANS export')
    hist=frames['rans_drift_series.csv']; surface=frames['rans_surface_forces.csv']; integrals=frames['rans_transfer_integrals.csv']
    summaries=[]; pair_summaries=[]; events=[]; turb_ends=[]
    for model in ['SA','SST']:
        all_steps=sorted(hist.loc[hist.model.eq(model),'step'].unique())
        wide=pd.DataFrame({'step':all_steps}); wide['tconv']=wide.step*1e-4
        turb=frames[f'rans_turb_state_norms_{model}.csv']
        twide=pd.DataFrame({'step':sorted(turb.loc[turb['case'].str.contains('_Drift_'),'step'].unique())})
        twide['tconv']=twide.step*1e-4
        for case, h in hist[hist.model.eq(model)].groupby('case',sort=True):
            _, method, seed=case_info(case); ident=f'{method}_s{seed}'
            h=h.sort_values('step'); s=surface[surface['case'].eq(case)].sort_values('step')
            if not np.array_equal(s.step,h.step): raise ValueError('Surface/history times differ')
            for coeff in ['dCL','dCD','dCMz']:
                wide[f'{ident}_{coeff}']=h[coeff].to_numpy()
            for coeff in ['dCDp','dCDv']:
                wide[f'{ident}_{coeff}']=s[coeff].to_numpy()
            t=turb[turb['case'].eq(case)].sort_values('step')
            twide[f'{ident}_eddy_rms_pct']=100*t.Eddy_Viscosity_band_rms_rel.to_numpy()
            cols=[c for c in t if c.endswith('_band_rms_rel')]
            end=dict(model=model,method=method,seed=seed,step=int(t.step.iloc[-1]),n_band=int(t.n_fixed_wallband.iloc[-1]))
            for col in cols: end[col]=float(t[col].iloc[-1])
            turb_ends.append(end)
            row=dict(model=model,method=method,seed=seed,n=len(h),first=int(h.step.min()),last=int(h.step.max()))
            for coeff, data in [('CL',h),('CD',h),('CMz',h),('CDp',s),('CDv',s)]:
                stat=ols(data.step.to_numpy(),data['d'+coeff].to_numpy())
                for name,v in stat.items(): row[f'{coeff}_{name}']=v
            summaries.append(row)
            I=integrals[integrals['case'].eq(case)]
            state_steps=sorted(I[I.which.eq('state')].step.unique())
            visible=[k for k in state_steps if k+1<=h.step.max()]
            events.append(dict(case=case,exported_history_pairs=len(state_steps),pairs_with_poststep_in_window=len(visible),
                               excluded_terminal_state_steps=[int(k) for k in state_steps if k not in visible]))
        write_csv(gen/f'rans_drift_{model}.csv',wide)
        write_csv(gen/f'rans_turbulence_{model}.csv',twide)
        for seed in [1,2]:
            n=wide[f'NN2_s{seed}_dCD'].to_numpy(); g=wide[f'GP2_s{seed}_dCD'].to_numpy()
            pair_summaries.append(dict(model=model,seed=seed,NN2_rms=np.sqrt(np.mean(n*n)),GP2_rms=np.sqrt(np.mean(g*g)),
                                       paired_rms=np.sqrt(np.mean((n-g)**2)),paired_max=np.max(np.abs(n-g))))
    summary=pd.DataFrame(summaries).sort_values(['model','method','seed'])
    # Display-scaled columns keep the manuscript tables compact. All unscaled values remain above.
    for field in ['CL','CD','CMz','CDp','CDv']:
        summary[f'{field}_trend_micro']=summary[f'{field}_trend']*1e6
        summary[f'{field}_se_micro']=summary[f'{field}_se']*1e6
    write_csv(gen/'rans_drift_summary.csv',summary)
    write_csv(gen/'rans_pair_summary.csv',pd.DataFrame(pair_summaries))
    write_csv(gen/'rans_turbulence_endpoint.csv',pd.DataFrame(turb_ends))
    integral_rows=[]
    valid=[]
    for item in events:
        I=integrals[integrals['case'].eq(item['case'])].copy()
        # Export names identify the latest level. The older level is one index earlier.
        I['event_step']=I.step+np.where(I.which.eq('state'),1,2)
        max_step=int(hist[hist['case'].eq(item['case'])].step.max())
        I=I[I.event_step<=max_step]
        valid.append(I)
    valid=pd.concat(valid,ignore_index=True)
    for model, fields in [('SA',['Density','Energy','Nu_Tilde']),('SST',['Density','Energy','Turb_Kin_Energy','Omega'])]:
        for field in fields:
            row={'model':model,'field':field}
            for method, token in [('NN2','_NN_'),('GP2','_ConsGalerkinProj_')]:
                a=valid[valid.model.eq(model)&valid.field.eq(field)&valid['case'].str.contains(token,regex=False)]
                y=a['rel_defect_(tgt-src)/src'].abs()
                if len(y)!=156: raise ValueError(f'Expected 39 events x two histories x two seeds: {model} {field} {method}: {len(y)}')
                row[method+'_max']=y.max(); row[method+'_median']=y.median(); row[method+'_n']=len(y)
            row['max_ratio_NN2_GP2']=row['NN2_max']/row['GP2_max']
            integral_rows.append(row)
    integral_summary=pd.DataFrame(integral_rows)
    integral_summary['quantity']=integral_summary.field.map({'Density':r'$\rho$', 'Energy':r'$\rho E$', 'Nu_Tilde':r'$\rho\widetilde\nu$', 'Turb_Kin_Energy':r'$\rho k$', 'Omega':r'$\rho\omega$'})
    write_csv(gen/'rans_integral_summary.csv',integral_summary)
    # Report-derived setup and limiter metadata are kept distinct from measured CSVs.
    report=(root/'MDFiles'/'NACA0012_STUDY_REPORT.md').read_text()
    clip=read_markdown_table(report,'| case | setting | field | transfers | floor |')
    clip=clip[clip['case'].str.contains('Drift') & clip.field.eq('Turb_Kin_Energy')].copy()
    clipping=[]
    for _,r in clip.iterrows():
        clipping.append(dict(model='SST',method='GP2' if 'ConsGalerkinProj' in r['case'] else 'NN2',
             seed=2 if 'Seed2' in r['case'] else 1,exported_history_transfers=int(r.transfers),with_raises=int(r['with raises']),
             max_nodes=int(r['max nodes']),min_before=first_number(r['lowest before']),floor=first_number(r.floor)))
    write_csv(gen/'rans_clipping_report.csv',pd.DataFrame(clipping))
    # Native comparisons: their deltas are rounded source summaries, not complete time series.
    native=[]
    for _,r in frames['rans_deformed_reference.csv'].query("what == 'steady'").iterrows():
        native.append({c:first_number(r[c]) if c!='model' else r[c] for c in ['model','CL_original','CD_original','CDp_original','CDv_original','dCL','dCD','dCMz']})
    write_csv(gen/'rans_native_summary.csv',pd.DataFrame(native))
    # One-shot views use actual recorded samples only; shorter windows end in NaN, not extrapolation.
    one=frames['rans_oneshot_series.csv']; lags=frames['rans_oneshot_lags.csv']
    onesummary=[]
    for model in ['SA','SST']:
        wide=pd.DataFrame({'lag':np.arange(int(one.steps_after_event.max())+1)})
        for case, a in one[one.model.eq(model)].groupby('case'):
            meta=lags[lags['case'].eq(case)&lags.coeff.eq('CL')].iloc[0]
            tag=('d005' if meta.d_in==.05 else 'd010' if meta.d_in==.1 else 'd020')
            if '_SF2_' in case: tag+='_SF2'
            if '_Long_' in case: tag+='_long'
            a=a.sort_values('steps_after_event')
            series=a.set_index('steps_after_event').dCL
            wide[tag]=wide.lag.map(series)
            row=dict(model=model,band=tag,d=meta.d_in,event=int(meta.event),window=int(meta.window),
                     pre_floor=meta.pre_floor,peak=meta.peak,peak_lag=int(meta.peak_lag),
                     lag50=meta.lag_50pct,predicted=meta.predicted_lag,censored=meta.censored)
            onesummary.append(row)
        write_csv(gen/f'rans_oneshot_{model}.csv',wide)
    one_table=pd.DataFrame(onesummary)
    one_table['case_label']=one_table.band.map({'d005':'0.05','d010':'0.10','d020':'0.20','d010_SF2':'0.10 / SF2','d020_long':'0.20 / long'})
    one_table['censor_label']=one_table.censored.map({'yes':'C','no':'--'})
    write_csv(gen/'rans_oneshot_summary.csv',one_table)
    probe=frames['rans_pressure_probe.csv'].copy()
    ulp=2**-24
    probe['peak_in_ulp']=probe.peak/ulp;probe['below_three_ulp']=probe.peak<3*ulp
    write_csv(gen/'rans_pressure_probe_audit.csv',probe)
    # Compact surface summary at the end of the drift window.
    ss=[]
    for case,s in frames['rans_surface_norms.csv'].groupby('case'):
        if '_Drift_' not in case: continue
        model,method,seed=case_info(case); r=s.loc[s.step.idxmax()]
        ss.append(dict(model=model,method=method,seed=seed,step=int(r.step),rms_dCp=r.rms_dCp,max_dCp=r.max_dCp,
                       rms_dCfx=r.rms_dCfx,max_dCfx=r.max_dCfx,maxYplus_ref=r.maxYplus_ref,maxYplus_case=r.maxYplus_case))
    write_csv(gen/'rans_surface_endpoint.csv',pd.DataFrame(ss))
    # Dynamic manuscript numerals: ensure text follows changed input CSVs as well as plots.
    def cmd(name: str, value: float, fmt: str='.2f') -> str:
        return '\\newcommand{\\'+name+'}{'+format(value,fmt)+'}\n'
    macros='% Generated from RANS_PostProcessing CSVs; do not edit.\n'
    macros+=cmd('RansCdRmsMin',summary.CD_rms.min()*1e5)+cmd('RansCdRmsMax',summary.CD_rms.max()*1e5)
    macros+=cmd('RansPairRmsMin',pd.DataFrame(pair_summaries).paired_rms.min()*1e6)+cmd('RansPairRmsMax',pd.DataFrame(pair_summaries).paired_rms.max()*1e6)
    macros+=cmd('RansViscTrendMin',abs(summary.CDv_trend).min()*1e7)+cmd('RansViscTrendMax',abs(summary.CDv_trend).max()*1e7)
    macros+=cmd('RansEddyMin',pd.DataFrame(turb_ends).Eddy_Viscosity_band_rms_rel.min()*100)+cmd('RansEddyMax',pd.DataFrame(turb_ends).Eddy_Viscosity_band_rms_rel.max()*100)
    macros+=cmd('RansUnderResolvedShells',int(probe.below_three_ulp.sum()),'d')+cmd('RansProbeShells',len(probe),'d')
    macros+=cmd('RansMaxYplus',frames['rans_surface_norms.csv'].maxYplus_case.max())
    write_text(gen/'rans_numbers.tex',macros)
    digest=json.loads((audit_dir/'rans_input_sha256.json').read_text())
    changed=[name for name,h in digest.items() if not (root/name).exists() or hashlib.sha256((root/name).read_bytes()).hexdigest()!=h]
    # Intentional future input updates are allowed, but must be disclosed against this delivered snapshot.
    audit=dict(source_hashes=hashes,source_tables=audits,event_accounting=events,
               duplicate_policy='Remove exact duplicate rows in derived views; fail on conflicting physical keys.',
               integral_policy='Use exported relative-defect field, not subtraction of rounded inventory columns; both histories, 39 events, two seeds.',
               statistics='OLS on unique physical steps; standard errors are unadjusted and descriptive, not inferential.',
               probe_count=len(probe),probe_below_three_ulp=int(probe.below_three_ulp.sum()),
               pressure_ulp=ulp,short_run_convective_duration=.04,
               old_and_raw_csv_snapshot_count=len(digest),inputs_differing_from_delivered_snapshot=changed,
               bibliography='References.bib',raw_simulations_rerun=False)
    write_text(audit_dir/'rans_data_audit.json',json.dumps(audit,indent=2,allow_nan=False)+'\n')
    print(f'Built RANS views; preserved {len(digest)} CSV inputs; {len(changed)} differ from archived snapshot.')
    print(f'Removed only exact duplicates in views: '+str({k:v['exact_duplicates_removed_in_view'] for k,v in audits.items() if v['exact_duplicates_removed_in_view']}))
    print(f'{probe.below_three_ulp.sum()} / {len(probe)} pressure-shell peaks lie below three source-field ULPs.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args()
    build(args.root.resolve())
