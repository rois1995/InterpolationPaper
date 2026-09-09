# Scientific audit and substantive revision

## Editorial decision

**The study has a coherent paper, but the supplied evidence does not yet support a paper claiming a completed assessment of unsteady RANS remeshing.** It supports a controlled Euler investigation of transfer accuracy, accumulated remeshing effects, BDF2 history forcing, and aerodynamic response. The most distinctive evidence is the mixed-history NACA experiment, supplemented by the paired raw/global-correction comparisons and the pressure-propagation diagnostics. The principal contribution is not invention of the interpolation operators or the BDF2 coefficients.

The revised manuscript preserves the intended RANS motivation but distinguishes it from the completed Euler evidence. Its title is **Mesh-to-Mesh Transfer Errors in BDF2 Unsteady Flow Simulations**. This is not a change to the long-term scientific objective: it prevents the current draft from claiming a turbulence-model assessment that has not been performed. A genuinely unsteady RANS case should become the final application before submission under the intended scope.

The rewrite reduces the supplied 42-page manuscript to 28 pages in the same 11-point, one-inch-margin format. It does not obtain that reduction by shrinking the body font or deleting contrary results. The main argument now occupies approximately 16 pages; secondary tables, convergence checks, representation analysis, and verification details follow in appendices. Page counts are layout measures, not estimates of journal typeset length.

## Basis of the audit

The reviewed project is the supplied `interpolation_paper_naca_results_package.zip`, including `paper.tex`, `References.bib`, 39 existing input/generated CSVs, the processing scripts, and the Riemann and NACA reports. The reusable scientific-review prompt was read in full. The reports were treated as evidence summaries, not as independent verification of their own explanations.

The numerical checks were performed on the available tabulated values. The actual flow solver, mesh files, raw restart fields, residual vectors, complete histories, and case-directory code copies are not present. I did not rerun simulations, inspect their implementation, or independently reproduce raw-field residual, pressure-probe, or conservation measurements. In particular, the named CSVs described by the reports are not necessarily the same as the report-transcribed CSVs supplied with the manuscript.

External primary literature was used to verify attribution and to assess the paper's novelty and journal fit. All additions are identified below. Source reports are retained unchanged; scientific corrections are made explicitly in the manuscript and documented here, rather than silently rewriting the evidence record.

## 1. Scientific corrections implemented

### 1.1 A native-target reference is not an exact interpolation-error reference

The previous draft repeatedly stated that comparison against a native target trajectory isolates the interpolation error. That is too strong. For a transferred source state, the discrepancy

`epsilon = I_AB(U_A) - U_B`

still contains the native source/target discretization mismatch. Using the target trajectory is essential for a common evolution operator, but does not remove that mismatch. The revised paper calls this the **target-referenced restart discrepancy**. A source/no-remeshing reference in the repeated-remeshing experiments answers a different operational question: the combined effect of changing discretization and transferring states.

The error decomposition has also been corrected for nonlinear transfer. Limiting and bounded redistribution generally prevent writing `I(Pu + e) = I(Pu) + I(e)`. The new exact decomposition uses differences of the complete nonlinear map and reduces to the familiar additive expression only when the map is linear.

**Consequence:** the higher-order NACA norm cluster is compatible with a common native mismatch, but its existence does not quantify that mismatch or prove it dominates. A unique GP3 accuracy ranking is not supported; BAR2 has the smallest tabulated target-referenced norm.

### 1.2 The BDF2 derivation has a defined domain of validity

The derivation now applies to a fixed target mesh between discrete remeshing events. A general ALE calculation with moving control-volume measures cannot be represented simply by placing every geometric effect inside a spatial residual while retaining a constant mass matrix in the temporal term. Moving-volume histories and GCL consistency need their own formulation.

The exact common-trial history-residual identity is distinguished from the local state-response linearization. A locally Lipschitz residual Jacobian, a locally invertible implicit step, small perturbations, and sufficiently small inner-solve residuals are required for the quadratic remainder statement. Merely saying that a residual is differentiable is insufficient for an `O(||delta||^2)` remainder.

Remaining algebraic errors now appear explicitly as `dt*(r_hat-r)` in the perturbation equation. This prevents solver error from being silently assigned to interpolation or to the nonlinear remainder.

### 1.3 A force-slope discrepancy is not the first-step discrepancy divided by dt

The correct discrete slope difference is

`(Delta G[n+1] - Delta G[n])/dt`.

The previous expression omitted the discrepancy already present at the transferred latest history. The revised expression includes both the response at the new step and the output gradient acting on `epsilon[n]`. Therefore, an unconditional first-order-transfer/finite-kink argument based only on `Delta G[n+1]/dt` has been removed. Cancellation can matter, including in the steady full-history limit.

### 1.4 Galerkin projection is not automatically cellwise conservative or spatially local

For an unconstrained continuous target test space containing constants and coordinates, exact Galerkin orthogonality preserves the consistent zeroth and first moments; P2 additionally contains quadratic coordinate functions. It does not impose a conservation equation on every finite-volume cell, whose indicator is generally outside the continuous space.

The consistent mass matrix is sparse, but its inverse is not generally local. A projection is the identity when the complete input field is already in the complete target space; identical source and target geometry on a subregion does not establish an identity transfer there. Global limiter redistribution adds another nonlocal operation. The statement that an annular mesh change confines the transfer error exactly to the annulus has therefore been removed.

The Riemann conclusion is narrowed to the supported result: **enforcing a total integral does not determine the spatial transport error**. It does not prove that only a locally conservative method can maintain the transition position. This is particularly important because higher-order NN/BAR variants also perform well in the frequency comparison.

### 1.5 Vertex-only P2 output does not have a universal third-order lumped defect

The source material associates discarded P2 edge information with an `O(h^3)` lumped-integral discrepancy. This is not generally true. The revised appendix supplies an explicit analytical counterexample: the exact P2 representation of `f(x)=x^2` on `[0,1]`, sampled at uniform vertices and integrated with P1 nodal weights, gives

`Q_h(f) = 1/3 + h^2/6`.

The nodal values and projected field are exact, but their vertex-lumped integral has a second-order error. The same counterexample extends to `f(x,y)=x^2` on a uniformly triangulated unit square. Source lifting, target restriction, and weight compatibility must be distinguished. This analytical check does not change any experimental CSV.

**Consequence:** the order of a global correction is not automatically the order of the underlying point interpolant. Preserving pointwise order requires a correction small enough in that norm, and a small total defect does not itself bound the maximum local change when admissible capacity becomes small. Correcting the stored nodal field also does not preserve all properties of a detached, unmodified full P2 field.

### 1.6 Limiter claims now match the actual formulas

BAR3 has a half-gradient correction, not an explicit Taylor Hessian term. Setting its correction factor to zero recovers BAR2. For NN3, damping both Taylor terms recovers NN1. A bound derived from linear increments does not automatically bound a quadratic at all evaluation points.

A raw L2 projection's best-approximation and moment properties cannot be assigned unchanged to its clipped/redistributed output. Source-vertex bounds need not contain every smooth quadratic subcell extremum. Componentwise bounds on density, momentum, and total energy do not alone imply positive pressure. The revised manuscript treats these as properties to verify on the complete returned state rather than unconditional guarantees.

## 2. Results audit

### Recomputed quantities

| Check | Recalculated from available values | Editorial consequence |
|---|---:|---|
| NN2/GP2 two-seed mean numerical-drag trend | 5.8333 | Report approximately 5.8, not a universal factor of six. |
| NN1 GC change in full-history lift response | +0.08864% | Remove “unchanged in the fourth significant figure.” |
| NN1 GC change in state discrepancy norm | −0.01085% | Small, not exactly zero. |
| NN1 GC change in first-step state norm | −0.00890% | Small, not exactly zero. |
| BAR2 GC change in lift response | −0.002737% | Supports weak sensitivity for this pair. |
| Higher-order NACA discrepancy-norm spread | 2.8835% | Does not support a unique GP3 ranking. |
| NN1 peak/first-step lift discrepancy | 4.6782 | First-step response does not bound the later peak. |
| Mean pressure speed over six thresholds | 1.18367 | Compatible with the acoustic scale. |
| Sample SD across those threshold speeds | 0.07916 | Detector sensitivity, not an independent-sample confidence interval. |
| Pressure speed range | 1.0905–1.2788 | Main paper reports the range and the nature of the spread. |
| Riemann fitted orders | 1.0101, 2.4001, 2.9873, 1.9734, 3.2026, 2.1217, 2.7915 | Consistent with the rounded source table. |

The static generated CSVs were rebuilt into a separate temporary directory; all nine regenerated CSVs matched the supplied generated files byte-for-byte. All 39 pre-existing CSVs in the working project match their original SHA-256 hashes.

### 2.1 NACA residual diagnostics contain an invalid supporting statistic

The reported normalized cosines exceed one, in the range approximately 1.0000008–1.0000017. A correctly normalized cosine cannot exceed one. The underlying residual comparison may still be sound, but this statistic cannot be displayed as evidence of extraordinary agreement.

The paper retains the reported relative residual-discrepancy ranges and states their export limitation. Cosine columns are preserved in the original CSV but omitted from the printed table. They must be recomputed from the original vectors, in a consistent metric and preferably double precision. Do not simply clamp them to one: that would hide the diagnostic defect.

### 2.2 Native NACA relaxation is not random noise

The native target trajectory changes smoothly after restart despite a deeply converged steady density residual. Its 51-step coefficient span is not a stochastic noise floor or an independently established algebraic-error bound. Subtracting the native trajectory is appropriate, but does not demonstrate that every perturbed physical step is converged below the observed signal.

The older-history-only result remains strong: the latest history is identical to the native one, whereas the older history changes and the response is nonzero with the expected sign. The rewrite emphasizes this controlled contrast rather than a “10–13 times the noise floor” significance claim.

The measured norm ratios support approximate directional linearity, but vector superposition defects remain 1.8–3.5%. The statement that no amplitude study is needed to establish a general linear regime has been removed. There is no direct evaluation of the target Jacobian/adjoint prediction. With only one small time step, agreement with the history ratios can also be dominated by the temporal term; a time-step sensitivity check would test more of the spatial resolvent than those ratios alone.

### 2.3 Global correction preserves the source total, not the native target total

The NACA integral-difference column is target-referenced. Its corrected common value of −8.12114e−4 represents the native source/target integral difference; it is not failure of a correction intended to preserve the source total. The revised text and table caption identify this reference explicitly.

The small response changes are useful evidence, but only for the supplied pairs and settings. They do not establish that global correction never affects force response.

### 2.4 Repeated-load statistics do not establish zero-mean lift or moment error

GP2 lift slopes change sign across two seeds. Both NN2 lift slopes are positive and different in magnitude. NN2 moment slopes are closely matching and negative. These data contradict the blanket explanation that lift and moment errors cancel while drag alone accumulates.

The drag trend is a reproducible two-realization observation, not an ensemble distribution. Reported ordinary regression standard errors are retained as such, not as autocorrelation-robust confidence intervals. No entropy budget or controlled separation of mesh-operator changes establishes dissipation as the unique mechanism. The manuscript no longer claims that exact conservation eliminates all other causes.

### 2.5 Acoustic interpretation is supported, but its precision was overstated

Six thresholds are applied to one pressure event, not to six independent events. Their spread measures sensitivity to the chosen detector. The three force-onset radii give negative fitted intercepts of roughly 8.5–18 steps in magnitude, significant relative to the earliest onsets. High R² from three points does not prove a zero-intercept law or uniquely distinguish wave transport from all alternatives.

A fraction-of-own-peak detector is invariant to pure amplitude scaling. Thus the amplitude-control result is a useful extraction-stability check but cannot by itself exclude a threshold explanation. The pressure speeds and the separation from the much slower convective scale are the stronger evidence for acoustic transport. Finite-Mach directional characteristic speeds need not equal a single freestream sound speed exactly.

### 2.6 Riemann frequency/amplitude results remain limited by tracker provenance

The maximum-gradient and 50%-crossing diagnostics are separated. More importantly, the report records substantially different NN1 crossing estimates under different samplings for the same nominal configuration. This can compromise an amplitude factor if compared extractions do not use identical settings.

The requested frequency/amplitude comparison remains in the main paper, with the uncertainty attached to that figure and discussion. The strongest conclusion remains the paired raw-versus-GC comparison under common sampling, not “280 times better,” “position immunity,” or a universal tripling law. The full tables remain in the appendix.

Positive accumulated density/energy values are described as gains, not losses. The GP3 width changes cited in the report's interpretation come from a frequency comparison; there is no matching GP3 amplitude pair in the supplied table. No GP3 amplitude conclusion is retained.

### 2.7 Moment/limiter measurements are supporting observations, not a mechanism proof

A domain-density centroid is not a particular transition crossing. The 10−7 to 10−5 lumped moment changes are not machine precision. The GP2 global excursion being below a global tolerance does not prove that every local bound is inactive. No resolved per-cycle or accumulated limiter-on/off width budget is supplied. The paper therefore uses these observations to motivate an interpretation, without claiming a causal decomposition or displacement scaling law.

### 2.8 Summary files do not support continuous curves or every intermediate-time claim

`data/naca_force_history.csv` contains first, maximum, and final values, not a time series. `data/naca_surface_delta_cp.csv` contains maxima and their locations, not complete surface distributions. The rewritten paper uses tables or isolated measured values accordingly. It does not draw smooth histories or Cp curves through quantities that were not supplied.

## 3. Structural and writing changes

The narrative is now: why remeshing matters for unsteady RANS; what each transfer actually computes; how transferred histories enter a fixed-target BDF2 problem; how the references and experiments are constructed; and what the static, Riemann, and NACA evidence demonstrates. Mathematical conditions are defined once and referenced in the results.

The manuscript removes repeated statements of the paper's intentions, multiple lists of nominal properties, duplicated force/adjoint interpretations, planned but unperformed measurements described as completed, and broad “local conservation” or “immunity” conclusions. It retains one main smooth convergence figure, one discontinuous convergence figure, the representative Riemann stress comparison, the NACA mixed-history response, response ratios, repeated drag, and delayed force onset. Conservation curves, repeated-transfer details, full first-step matrices, moment and limiter audits, and projection patch tests are in appendices.

The BDF2 coefficients are treated as established algebra used to create a diagnostic experiment. They are not promoted as an independently novel time-integration theorem. The substantive contribution must be the controlled evidence and its consequences for remeshed CFD.

No obsolete conservative interpolation family is reintroduced. BAR3 is described only as the validated formulation used for the reported result. Source reports are not rewritten to conceal or reconcile their differing interpretations.

## 4. Required before submission under the intended RANS scope

### Priority 1 — Add the unsteady RANS evidence that the stated objective requires

One genuinely unsteady RANS application is more important than another static interpolation table. It should include the complete transported turbulence state and both BDF2 histories; pressure and viscous force contributions; identical forcing and comparison times; and controlled spatial, temporal, and inner-iteration errors. Demonstrate the effect on an actual engineering output such as load amplitude, phase, mean force, or a separation-sensitive quantity. Retain the Euler cases as mechanism-isolating verification rather than replacing them.

For causal attribution, separate the same mesh/time schedule with different transfer choices from comparisons that also change adaptation decisions. If adaptation follows each method's evolving solution, differing mesh schedules are an additional feedback effect and must be identified. At least a prescribed-mesh replay comparison is valuable.

### Priority 2 — Establish that the physical-step solve is not determining the measured response

Archive steady and dual-time configurations; export the exact residual definitions, signs, norms, precision, and complete temporal terms; demonstrate a native restart plateau or explain the native relaxation; and tighten physical-step tolerances until the target-referenced response stabilizes. Recompute the invalid cosine diagnostic from raw vectors. An exceptionally small steady density residual does not replace these checks.

### Priority 3 — Archive reproducible Riemann and NACA measurements

Supply the original case-generated CSVs, full reference and response histories, pressure-shell locations/onsets, volume weights, state and residual vectors needed for the checks, and versioned extraction scripts. Specify the Riemann initial states and physical boundary conditions rather than marker names alone. Resolve the conflicting NN1 tracker sampling. Where full raw fields cannot be distributed, archive sufficient reduced data and exact reduction code.

### Priority 4 — Close the linear-response and representation questions

Use a scaled perturbation sequence to determine whether superposition/model errors decrease as expected, with an independently established solver floor. Add non-collinear history perturbations from an unsteady trajectory. A useful inexpensive algebra check is a nonzero cancellation pair `epsilon[n]=epsilon/4`, `epsilon[n-1]=epsilon`, for which `eta=0`; in the stated fixed-mesh formulation the same unique solved step should result when all other inputs are identical. This is a proposed check, not a reported result.

For GP, separately export source-lift and target-full consistent integrals, final vertex-lumped totals, and pre/post-limiter changes. Check exact polynomial reproduction on irregular meshes and distinguish source lifting from target restriction. Record whether active limiting or boundary enforcement changes admissible test functions.

### Priority 5 — Make an engineering method recommendation supportable

Measure total transfer wall time, geometry/reconstruction/solve costs, memory, and extra physical-step iterations. Use enough mesh realizations to assess the actual variability of the claimed effects; two seeds are a pilot, not a universal sufficiency threshold. Cost and robustness could change the practical ranking even when accuracy favors a projection.

Additional three-dimensional testing and re-discretized surfaces would broaden the paper, but are not substitutes for the missing RANS application and solver-error controls. The most efficient next work is to resolve the central claim, not to expand every auxiliary parameter sweep.

## 5. Literature and attribution

The original 40 bibliography entries and their keys are preserved. Four verified records are appended to `References.bib`:

- `MaddisonHiester2017`: Optimal constrained interpolation; DOI **10.1137/15M102054X**. Shows why integral preservation and optimal L2 approximation do not preserve every discrete property.
- `MaddisonCotterFarrell2011`: Balance-preserving mesh transfer; DOI **10.1016/j.ocemod.2010.12.007**. Relevant precedent for preserving a discrete steady balance rather than only field accuracy.
- `MaddisonCotterFarrell2013`: Corrigendum; DOI **10.1016/j.ocemod.2013.04.007**. Included with the original record rather than ignoring its corrections.
- `SauvageAlauzetDervieux2024`: Space/time mesh adaptation for implicit transient RANS; DOI **10.1016/j.jcp.2024.113389**, JCP 519, 113389. Places the intended application in the relevant adaptation literature.

Metadata and the cited scope of these works were checked using publisher pages or authors' institutional repositories. This was a targeted attribution/novelty check, not an exhaustive literature review or a verification of every assertion in every cited publication. The retained bibliography may contain uncited records; it remains the user's master file rather than being destructively pruned.

## 6. Reproduction of the revision

`make` builds the PDF from `paper.tex`, `References.bib`, and the existing CSV paths. `make audit` recomputes the numerical checks. The preserved CSV hashes are stored in `audit/original_csv_sha256.json`; `python scripts/audit_results.py --check-preservation` verifies this reviewed snapshot. Intentional future data updates will naturally require a new provenance baseline rather than alteration of this one.

The delivered audit JSON gives the actual calculations, and the build-validation JSON records successful compilation, citation/reference checks, and source/data preservation. All pages were rendered and inspected, including the revised equations and figures. A misleading baseline in an inherited logarithmic conservation bar chart was replaced by a point plot without changing its data. No source or font assets outside the project are needed other than the normal installed LaTeX packages.

The revision is substantially more coherent and scientifically defensible. **It is not a certification that the numerical campaign is complete or submission-ready.** The missing RANS evidence and raw-data/solver checks are research obligations; they cannot be discharged by more polished prose.
