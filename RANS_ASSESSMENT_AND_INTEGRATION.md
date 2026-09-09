# Assessment and integration of the fixed-body NACA 0012 RANS results

## Editorial decision

The RANS results belong in the paper and complete a **fixed-body, time-accurate RANS application** of its transfer/remeshing question. They do not reproduce the strongest Euler method ranking or the Euler acoustic-delay measurement. That contrast is useful evidence: integral preservation, localized state changes, and aerodynamic response must be assessed separately. The application should not be framed as another demonstration that projection always outperforms pointwise transfer.

The main new result is **substantially smaller recorded turbulence-inventory defects with GP2, without a correspondingly clear aggregate load advantage over NN2 in the tested short RANS runs**. A second result is that the fixed near-wall vertices develop a measurable eddy-viscosity discrepancy despite the protected geometry. Neither observation establishes the isolated causal contribution of turbulence interpolation.

The manuscript now includes the RANS setup, a paired load-history figure, and an instantaneous transfer-integral table. Detailed slopes, floor enforcement, fixed-band turbulence histories, and unsuccessful propagation diagnostics are in an appendix. The abstract, introduction, evidence statement, and conclusions have been revised consistently. The mathematical formulation and the existing controlled Euler BDF2 evidence are retained; the new RANS calculations are not relabeled as mixed-history validation.

No moving-body case is required to support this scope. There is also no reason to repeat the entire existing Euler parameter matrix with RANS. Remaining checks are limited to the accuracy and provenance of the present fixed-body results.

## Sources and work actually performed

The base is the scientifically reviewed source project supplied as `interpolation_paper_scientific_review_package.zip`, not an earlier long draft. New inputs are:

- `MDFiles/NACA0012_STUDY_REPORT.md`, particularly Part III;
- `MDFiles/Scientific Manuscript Writing and Development Prompt.md`;
- all 14 CSVs and the accompanying reports/logs in `RANS_PostProcessing/`.

The source reports and CSVs are retained unchanged. I read the writing prompt, compared the report with the CSV schemas and values, checked physical-step uniqueness, recomputed slopes and paired load statistics, inspected integral/history indexing, screened the pressure-probe precision, and rebuilt the figures and tables from derived CSV views. The exact source hashes and automated audit are in `audit/rans_input_sha256.json` and `audit/rans_data_audit.json`.

I did **not** run the CFD solver, inspect original mesh/restart vectors, verify case source code, or reproduce the extraction of surface forces from VTU files: those assets and extractor scripts were not included. The new script reproduces the analysis of the reduced exports, not the simulation itself. Floor statistics are explicitly report-derived because the underlying nodal clipping logs were not supplied.

## 1. Scientific interpretation

### 1.1 Application scope and time horizon

The reported configuration is fixed NACA 0012 at Mach 0.1, incidence 6 degrees, and Reynolds number 6 million, with SA and SST on the same hybrid viscous mesh. It advances in physical time with BDF2 and replaces the mesh discretely. This is aligned with the paper's fixed-mesh-between-events derivation.

The RANS nondimensionalization has sound speed 1 and freestream speed 0.1. Consequently, 400 steps of size 0.001 span 0.4 code time units but only **0.04 chord-convective times**. The extended 800-step references span 0.08. The data are therefore evidence of short remeshing-induced transients and finite-window accumulated changes, not long-time RANS equilibrium statistics, long-time drift rates, or periodic-load accuracy.

The physical baseline is nominally steady, but the no-transfer reference is not exactly stationary after startup. Each case is compared against its own closure's reference at the same physical step. This removes the common reference trajectory from the plotted differences; it does not establish that every perturbed physical step is converged below the reported differences.

### 1.2 The large Euler GP2–NN2 load ranking does not carry over

Source: `rans_drift_series.csv`; analysis: `generated/rans_drift_summary.csv` and `generated/rans_pair_summary.csv`.

Total-drag RMS deviations from the native reference range from approximately 3.01e-5 to 3.43e-5. For each common random seed, the relative difference between the NN2 and GP2 RMS values is:

| Closure | Seed | 100 × (NN2 RMS / GP2 RMS − 1) |
|---|---:|---:|
| SA | 1 | −0.730% |
| SA | 2 | +0.950% |
| SST | 1 | −0.133% |
| SST | 2 | −0.347% |

There is no consistent ordering, and the largest difference is below 1%. This does **not** mean the individual trajectories coincide. The RMS of the step-aligned NN2–GP2 difference is 8.71e-6 to 8.84e-6; the maximum paired differences are about 3e-5. Equal aggregate RMS errors and equal time histories are different statements.

The revised main figure displays all four trajectories per closure. It does not replace them with an ensemble average, attach unsupported confidence intervals, or claim universal equivalence. In particular, the approximately 5.8 Euler numerical-drag trend ratio is not extended to RANS.

### 1.3 Pressure and friction drag must remain separate

Source: total history loads from `rans_drift_series.csv`; surface-integrated contributions from `rans_surface_forces.csv` after removal of exact duplicate snapshots.

The recomputed friction-drag slopes are negative for all eight runs, with magnitudes 2.86e-7 to 4.52e-7 per 100 steps. Larger total-drag oscillations are primarily pressure-related. There is no evidence here for the earlier broad explanation that interpolation necessarily produces a positive, sign-definite accumulated drag trend.

Surface integration is useful because six-decimal `forces_breakdown` values quantize small friction changes. It is **not** full-precision solver data: the underlying Cp and Cf fields are single precision and the integration rule differs from the solver's. The manuscript uses case-minus-reference differences formed with the same quadrature and does not identify more printed digits with a smaller physical uncertainty.

The drift data are strongly serially correlated: the lag-one correlation of the detrended total-drag residuals is approximately 0.79–0.84. The report's ordinary OLS standard errors and `resolved = |slope| > 2 SE` flag are not autocorrelation-robust significance tests. Recomputed unadjusted errors remain in the appendix as descriptive fit diagnostics. No inferential claim is built on them.

### 1.4 Better turbulence-inventory preservation does not guarantee a load improvement

Source: `rans_transfer_integrals.csv`; view: `generated/rans_integral_summary.csv`.

Pooling two seeds, both histories, and the 39 events followed by recorded physical evolution gives 156 samples per method/field:

| Closure | Inventory density | Maximum absolute relative defect, NN2 | GP2 |
|---|---|---:|---:|
| SA | density | 5.503e-11 | 2.191e-13 |
| SA | total energy density | 4.260e-11 | 2.871e-13 |
| SA | density-weighted SA variable | 2.484e-6 | 4.081e-9 |
| SST | density | 4.943e-11 | 2.270e-13 |
| SST | total energy density | 4.155e-11 | 3.030e-13 |
| SST | density-weighted k | 2.984e-6 | 2.593e-8 |
| SST | density-weighted omega | 8.146e-8 | 1.080e-9 |

These are instantaneous **source-to-target** changes, not native-target discrepancies and not cumulative changes caused by physical evolution. Turbulence production, destruction, and fluxes change turbulence inventories between transfers. Those changes cannot be counted as transfer-conservation errors.

The much smaller GP2 values are worth reporting, especially for the turbulence fields. They do not establish finite-volume cellwise conservation, uniform spatial accuracy, or a corresponding load advantage. The large farfield also makes density/energy domain totals insensitive to small localized perturbations. The manuscript keeps the source-relative integral test separate from the native-reference force test.

### 1.5 Protected coordinates do not imply an unchanged near-wall state

Sources: `rans_turb_state_norms_SA.csv`, `rans_turb_state_norms_SST.csv`.

The fixed inner set contains 12,233 vertices. At step 395 the eddy-viscosity discrepancy has nodal RMS equal to 0.38–0.44% of the reference maximum on this set, with similar values for NN2 and GP2. This is not a local relative error at every vertex, and it is not volume-weighted. These definitions are retained in the figure caption and methodology.

Small band-normalized omega errors must not be compared directly with k or eddy-viscosity errors: the normalization uses a very large near-wall omega maximum. The common-space comparison is valid because the vertices coincide, but it cannot identify whether the discrepancy comes from nonlocal transfer, the modified residual operator, physical propagation, or finite physical-step convergence.

The source's statement that the whole wall-layer quadrilateral region never moves does not follow from the stated d_w < 0.01 freeze cutoff, because that region reportedly extends to 0.053 chord. An additional topology-based mask could reconcile the statements, but no mesh/mask was supplied. The paper relies only on the explicitly identified fixed inner vertices and does not assert that all quadrilaterals are protected.

The actual maximum reported y+ is about 2.69. The manuscript does not describe the entire wall as y+ < 1 or infer near-wall grid convergence from the first spacing.

### 1.6 SST positivity handling is part of the returned method

Source: the clipping table in Part III.6 of the report, not an independent nodal clipping CSV.

The reported kinetic-energy floor is active in every exported drift history. Maximum raised-vertex counts are 423 and 415 for the GP2 seeds, and 19 and 13 for NN2. The lowest pre-floor values are −1.29e-7 and −1.62e-7 for GP2, and −2.38e-7 and −2.17e-7 for NN2. No omega-floor activation is reported. SA uses no analogous floor and has reported tiny negative density-weighted SA values of order 1e-16.

Counts alone do not rank error: locations, increments, and volume weights matter. These diagnostics demonstrate active nonlinear post-processing. They do not isolate its contribution to the inventory error or the wall stress. The appendix reports the facts without importing ideal projection guarantees into the clipped result.

### 1.7 Native mesh sensitivity is substantial but is not an instantaneous correction

Source: `rans_deformed_reference.csv`.

The first displaced mesh changes the independently converged lift by −1.891e-3 (SA) and −9.649e-4 (SST). The first ten-step GP2 continuations differ from the original reference by only +1.318e-5 and +1.150e-5 on average. Comparing equilibria and following a transferred state for ten steps are not the same experiment. The equilibrium difference must not be subtracted from an arbitrary early transient to manufacture a “pure interpolation” result.

This comparison directly supports the paper's distinction between native mesh discrepancy and target-referenced restart response. It does not prove that the native mismatch dominates every method difference. The new RANS archive does not contain a separate controlled mixed-history test or a target Jacobian prediction; those claims remain supported by the earlier Euler experiment only.

### 1.8 The new RANS records do not establish an acoustic travel-time law

Sources: `rans_oneshot_series.csv`, `rans_oneshot_lags.csv`, `rans_pressure_probe.csv`.

For the standard three launch bands, 50%-of-peak lift onsets are 3, 8, and 5 steps for both closures. The corresponding acoustic scales are 50, 100, and 200 steps. The first two cases have maxima near the end of the recording window; their peak-normalized thresholds are right-censored by that window. This is not a distance-proportional resolved arrival law.

Forty-four of 72 pressure-shell peak amplitudes lie below three source-field ULPs (pressure increment 2^-24). Many larger peaks cross at zero or one step. The three-shell fit implying a speed near 14 is consequently not a defensible acoustic measurement. It is excluded from physical interpretation, not silently transformed into agreement with the expected sound speed.

These data do not prove that acoustics are absent. They show that this detector/precision/initial-disturbance combination does not isolate acoustic propagation. A nonlocal projected or limited disturbance, the first implicit update, native operator change, and field precision are possible contributors. Their individual effects cannot be assigned from the reduced archive alone.

The doubled-amplitude cases have peak ratios of about 2.7 and nearly unchanged fraction-of-own-peak onsets. Such a detector is inherently insensitive to pure amplitude scaling; this observation cannot by itself prove propagation. The later-starting long case also has a different source state and event-dependent random mesh, so it is not treated as a controlled extension of the same signal.

## 2. Data corrections made only in derived analysis

### Exact duplicate snapshots

| Export | Raw rows | Exact duplicates | Unique rows |
|---|---:|---:|---:|
| rans_surface_forces.csv | 8,644 | 680 | 7,964 |
| rans_surface_norms.csv | 8,644 | 680 | 7,964 |
| rans_surface_profiles.csv | 40,700 | 13,320 | 27,380 |

All duplicated physical keys have identical values. The source surface-fit table counted 478 rows per drift case, although there are only 398 physical steps. Duplicates nonuniformly reweight the regressions. The new builder drops exact repeats, fails on any conflicting physical-key duplicate, and recalculates the fits. Original input bytes are preserved.

### Physical events versus history-file exports

Latest-level drift files are indexed 9, 19, ..., 399; older-level files 8, 18, ..., 398. The available force window ends at 399. The last pair has no subsequent physical-step response in that window. Main transfer statistics therefore pool 39 dynamically represented events, not 40 terminal/export pairs or 80 history files.

The report-only floor table still includes 80 history transfers, including the terminal pair; its caption states this explicitly. One-shot archives similarly contain an initial and terminal export pair but only one analyzed event. No numerical source value has been changed to conceal this distinction.

### Precision

The integral audit uses the provided high-resolution relative-defect column. Its very small values cannot be reproduced by subtracting the more coarsely printed large integrals. Force-history data, six-decimal force breakdown, and single-precision surface/volume data are explicitly distinguished. More digits in a derived quadrature do not imply more digits in its input field.

## 3. Placement within the manuscript

**Main text**

1. One compact fixed-body RANS setup subsection: closures, mesh, physical-time/normalization, transfer schedule, protected band, returned turbulence histories, and precision.
2. A four-panel load figure: SA/SST total drag and surface-integrated friction drag for both methods and seeds. Common axis ranges within each row and no smoothing.
3. A compact source-relative transfer-integral table. This is the strongest RANS-specific comparison and does not duplicate a load plot.
4. Short interpretation of near-wall turbulence changes and the absence of a resolved RANS acoustic law.

**Appendix**

- recalculated descriptive slopes and their unadjusted OLS errors;
- floor statistics, with report provenance;
- fixed-band eddy-viscosity histories;
- one-event peaks/onsets/censoring;
- duplicate and event accounting, field-precision limitations, and remaining metadata controls.

Full raw/derived tables, paired-statistic CSVs, and all available profiles remain in the source project rather than being printed exhaustively. The old static order figures and their reference slopes are unchanged. No RANS convergence plot is created because only one viscous mesh/time-step class is supplied.

## 4. What is established, provisional, or unavailable

| Component | Status | Consequence |
|---|---|---|
| RANS settings and model labels | Reported; configuration files not supplied | Use stated settings; archive exact model variants and inputs before submission. |
| Load histories and paired statistics | Established at the reduced-export level | Main application evidence, with finite-window/precision qualifications. |
| Smaller GP2 recorded turbulence-inventory defects | Established at the export level | Main quantitative result; not a universal load ranking. |
| Near-wall state differences | Established for the defined fixed-set metric | Supporting result; no causal decomposition. |
| Acoustic travel-time fit in RANS | Not established | Retain diagnostic limitation, do not report a physical speed. |
| Independent-target RANS mixed-history validation | Not performed in this archive | Do not extend Euler verification claims to an unperformed RANS experiment. |
| Physical-step and temporal convergence | Not demonstrated by iteration counts alone | Limits the resolution of small method differences. |
| Long-time/periodic/separated/moving-body conclusions | Outside supplied study | Not required for this paper's narrowed fixed-domain scope. |

## 5. Minimal remaining actions, not a new campaign

1. **Selected convergence control.** For one reference and representative perturbed runs, tighten the physical-step solve and refine dt while keeping remeshing at the same physical times. Establish which load and turbulence differences survive. This is the central numerical-confidence check; it is more valuable than adding another body motion or dozens of transfer variants.
2. **Archive the configuration evidence.** Supply exact solver/version, closure variant, freestream turbulence and thermal wall settings, displacement masks, and event indexing. Confirm whether any extra rule freezes all quadrilaterals. The fields indicate y+ up to about 2.69; do not substitute a nominal first spacing for an actual resolution check.
3. **Only if an additional RANS acoustic claim is desired:** export selected pre-transfer, post-transfer-before-update, and pressure time histories in double precision, with a sufficiently long common event window. This is optional for the present manuscript because the failed RANS arrival diagnostic is already interpreted honestly. It is not a reason to launch moving-body simulations.
4. **Only if recommending a production method:** add measured transfer cost and extra physical-step effort. Without them, report accuracy/response observations rather than an engineering optimum.

The existing scope can be completed without a pitching/plunging application. The decisive scientific issue is the accuracy of the finite-time comparisons already made, not expansion to a new physical problem.

## 6. Reproducibility and bibliography

Run `make` from the project root. The original static builder and CSV paths remain intact. `scripts/build_rans_data.py` reads `RANS_PostProcessing/*.csv` and the source clipping table, writes only `generated/rans_*` views and an audit JSON, and updates the LaTeX numerals/plots/tables through Makefile dependencies. `make rans-audit` reruns this process. Narrative interpretations must still be reviewed when source data change.

The canonical filename is `References.bib`. All 44 pre-existing entries are retained. Two model references are appended: `SpalartAllmaras1994` and `Menter1994`. Their publication metadata were checked against the primary-reference records on the Turbulence Modeling Resource SA and SST pages. These references identify the model families; they do not establish the exact compiled solver variants used in the supplied runs.

The final build/visual checks and source-preservation results are recorded separately in `audit/rans_build_validation.json`. No claims of CFD re-execution or raw-field reconstruction are made.
