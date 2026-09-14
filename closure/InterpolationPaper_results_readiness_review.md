# Complete results re-assessment and submission-readiness review

## Executive decision

The manuscript is now scientifically coherent, but **I would not submit it yet if the RANS aerodynamic comparison remains part of the headline contribution**. The Euler/Riemann/BDF2 evidence is strong enough for a substantive numerical-methods paper. The fixed-body RANS extension is valuable and should remain, but its load-level NN2-versus-GP2 comparison is not yet solution-verified: a large fraction of remeshed physical steps reaches the dual-time iteration cap, especially for SST. The current data therefore support RANS conclusions about the *returned transfer state and monitored integral defects*, but not a quantitative RANS ranking based on aerodynamic loads.

The most efficient path to submission is not another geometry, a moving body, or a larger interpolation matrix. It is a **small numerical-control closure campaign** on the existing RANS case, followed by one final manuscript pass.

## What was discarded and rebuilt

Per the requested editorial rule, all previous numerical interpretations were discarded except the smooth interpolation-order verification. The new Results and Conclusions were written from the recalculated repository CSVs rather than edited incrementally from the previous prose.

The new Results section is organized by scientific question:

1. smooth-field order verification (retained);
2. separation of boundedness, approximation order, and integral preservation;
3. repeated-remeshing spatial transport in the Riemann case;
4. two-history BDF2 restart response;
5. delayed and accumulated Euler aerodynamic response;
6. fixed-body RANS transfer and receiving-solver sensitivity.

The old Riemann position-drift narrative has been removed. The recomputed data show that the main stress-test error is **transition broadening**, not large feature displacement. The previous absolute values of 7--26 cells of *position* drift are no longer used.

The old RANS statement that 400 steps cover only 0.04 chord flow-throughs has also been removed. Executed-run startup data show a nondimensional physical step of 0.347224 and therefore 13.9 chord flow-throughs over 400 steps.

## Scientific assessment by evidence block

### 1. Smooth interpolation-order verification -- strong

This remains one of the cleanest parts of the paper. The test is independent of the flow solver and directly verifies the implementation of the nominal transfer orders. It is appropriately framed as consistency with first-, second-, and third-order smooth-field behavior rather than proof of universal asymptotic constants.

**Status:** submission-ready.

### 2. Discontinuity, boundedness, and global conservation -- strong supporting evidence

The recomputed discontinuous-field data show that the three commonly conflated properties are independent:

- NN2 unlimited: overshoot 0.8031, undershoot 0.4157, integral defect -5.26e-5;
- BAR3 unlimited: overshoot 0.1766, undershoot 0.1811;
- BJ limiting removes measured extrema but does not, by itself, restore the NN lumped integral;
- global correction reduces the lumped integral defect to approximately 1e-14;
- hard-limited GP2 is bounded and conservative to roundoff in the tested representation;
- vertex-only GP3 remains bounded but carries a representation-dependent lumped-integral defect.

The smooth one-transfer check also shows that the global correction changes L1 by less than 0.8% at the tested resolution while enforcing the monitored integral to roundoff. This is useful because it prevents the later failure of GC under repeated remeshing from being dismissed as a large one-step accuracy penalty.

**Status:** submission-ready as a verification/supporting subsection.

### 3. Repeated Riemann remeshing -- substantially stronger after recomputation

The recomputed common tracker gives a much more defensible result than the previous mixed-sampling analysis.

At 80 transfers, every method locates the feature within 0.14 cell. Under 400 transfers, the dominant degradation is profile width:

| condition | NN1 width change | NN1+GC | BAR2 | GP2 |
|---|---:|---:|---:|---:|
| 80 transfers, SF4 | 0.57 | -- | 0.25 | 0.00 |
| 400 transfers, SF4 | 21.41 | 18.36 | 3.24 | 0.00 |
| 400 transfers, SF2 | 26.36 | 25.65 | 22.91 | 0.06 |

At the same time the reported position errors remain sub-cell, including the aggressive cases. This changes the physical/numerical interpretation: low-order repeated transfer acts primarily as an accumulating diffusion mechanism rather than as a persistent advective displacement of the transition.

The raw-versus-GC comparison is especially strong. NN1+GC reduces the cumulative density-integral defect to approximately roundoff but leaves 18--26 cells of excess width. This is a clean counterexample to the idea that exact global conservation controls spatial remapping accuracy.

A key editorial correction is also necessary: these data **do not establish that only a locally conservative method can preserve the feature**. NN2, NN3, BAR3, GP2, and GP3 all remain much sharper than NN1/BAR2 in relevant comparisons. The robust conclusion is about the insufficiency of a scalar global constraint, not a unique mechanism assigned to Galerkin projection.

**Status:** strong and suitable for the main paper. One random mesh sequence remains a limitation but is not a blocker for the qualitative broadening contrast, which is very large.

### 4. Controlled NACA 0012 BDF2 history experiment -- strongest distinctive result

This is, in my view, the most publishable result in the paper.

The older-only configuration has an exact latest target history but a perturbed older history. It nevertheless produces a finite first-step state and load response. The latest/full state-norm ratio is 1.3330--1.3333 and the older/full magnitude is 0.3335--0.3341, matching the BDF2 coefficients. State-vector superposition errors are 1.8--3.5%; lift superposition errors are 0.28--1.49%.

This directly establishes a practical point that is easy to miss in remeshing workflows: transferring only, or diagnosing only, the latest state is insufficient for a multistep restart.

The new double-precision residual audit adds an important implementation-level result. NN/BAR variants satisfy the common-trial forcing identity at roughly 1e-7 relative error. GP2/GP3 have a larger whole-vector mismatch localized to wall-normal momentum at slip-wall nodes. The mismatch equals the wall-normal forcing introduced by projection, while the tangential/interior residual agrees at roughly 1e-7. This is not a failure of the BDF2 algebra; it shows that boundary admissibility is part of the actual transfer pipeline.

This should remain in the main text because it turns an apparent residual discrepancy into a physically meaningful implementation lesson.

**Status:** strong. No additional mixed-history simulation is required to support the core two-history conclusion.

### 5. Global correction in the NACA restart -- useful negative result

For NN1, the full-history lift response changes from 5.7873e-3 to 5.7924e-3 after global correction; BAR2 changes from -7.2852e-3 to -7.2850e-3. The response norm is similarly insensitive.

This complements the Riemann evidence: removing one scalar integral component does not necessarily remove the force-relevant spatial discrepancy.

**Status:** strong supporting result, but it should not be generalized to all flows or corrections.

### 6. Euler repeated-remeshing aerodynamic response -- useful but configuration-specific

Across two random realizations, the Euler GP2 drag trends are approximately 6.84e-6 and 6.38e-6 per 100 steps, versus 3.95e-5 and 3.84e-5 for NN2. The mean ratio is about 5.8.

This is a reproducible method difference in the tested Euler case. It is now explicitly presented as configuration-specific because the RANS cases do not reproduce the same ranking. The 400-step Euler window is also very short in convective units (about 0.047 chord flow-throughs), so this is an event-response experiment, not long-time aerodynamic statistics.

**Status:** useful in the main narrative, provided the limitation remains explicit.

### 7. Euler delayed propagation -- good supporting validation

The localized-transfer experiment remains useful for the paper's delayed-response argument. Pressure-shell fits on the long case give threshold-dependent speeds about 1.09--1.27, bracketing the nondimensional sound speed 1.183. These are multiple threshold extractions from a common event, so the spread is detector sensitivity rather than a confidence interval.

**Status:** credible supporting evidence. Keep concise; it is not the paper's principal contribution.

### 8. Fixed-body RANS transfer -- valuable, but the load ranking is not yet resolved

The executed-run provenance improves this section materially. The current setup is a fixed NACA 0012 at M=0.1, alpha=6 deg, Re=6e6 with SA and SST; 400 steps correspond to 13.9 chord flow-throughs, not 0.04.

The strongest RANS result is the transfer-event integral comparison. Maximum NN2/GP2 defect ratios are approximately:

- SA rho: 184;
- SA rho E: 164;
- SA rho nu_tilde: 247;
- SST rho: 177;
- SST rho E: 139;
- SST rho k: 37;
- SST rho omega: 37.

This is a clear demonstration that GP2 preserves the monitored transferred inventories far better than NN2 in the tested pipeline. For SST, the archived integral is evaluated before the primitive positivity floor, which must be stated.

The aerodynamic comparison is less decisive. Total-drag RMS deviations from the no-remeshing reference are about 6.1e-5--6.8e-5 for both methods. The NN2/GP2 ordering changes with closure and random realization; paired method differences have RMS about 0.8e-5--1.1e-5.

The problem is not that the load result is negative. The problem is that its numerical resolution is not established:

- SA repeated-remeshing cases reach the 800-inner-iteration cap in 144--171 of 398 physical steps;
- SST cases reach it in 277--333 of 398 steps;
- median terminal SST density residuals are only about log10(rms rho) = -11.1 to -11.3 against a requested -12 target;
- the repository contains an additional solver-control sensitivity ladder whose changes are comparable to, or larger than, the baseline remeshing load signal, but the L0/L1/L2 settings are not yet documented clearly enough in the manuscript to use it as a publication-quality convergence demonstration.

Therefore the current RANS data support the statement:

> GP2 produces much smaller monitored transfer-integral defects than NN2, while the present calculations do not resolve a corresponding aerodynamic advantage.

They do **not** support:

> NN2 and GP2 are aerodynamically equivalent, or one is more accurate than the other.

**Status:** valuable application evidence, but quantitative aerodynamic ranking is not submission-ready.

## Main strengths of the paper

### A. The paper has a genuine unifying question

The strongest conceptual contribution is not any individual interpolant. It is the decomposition of a remeshing event into distinct questions:

1. how accurately is the state represented on the target mesh?
2. are bounds/admissibility and monitored integrals preserved?
3. are all multistep histories transferred consistently?
4. what does the receiving discrete operator do with the resulting discrepancy?
5. does that response affect an aerodynamic output immediately or later?

The numerical campaign now maps cleanly onto that chain.

### B. The paper contains counterexamples rather than only favorable demonstrations

Several negative results materially improve credibility:

- exact global conservation does not repair repeated-remapping diffusion;
- high nominal order does not guarantee boundedness;
- the latest BDF2 history can be exact while the next state is still perturbed;
- the ideal Galerkin projection can generate a boundary-inadmissible momentum component that the solver then removes;
- a large advantage in transfer-integral preservation does not automatically become a resolved load advantage in RANS.

This is much stronger than a paper that merely ranks interpolation errors.

### C. The BDF2 experiment is unusually controlled

The full/latest/older decomposition is simple, interpretable, and directly tied to the mathematical formulation. It gives the manuscript a contribution beyond a generic remapping-method comparison.

### D. The recomputed Riemann analysis is now much cleaner

Using one common extraction removes the previous sampling-confounded position-drift argument. The new broadening result is both easier to defend and more relevant to repeated remapping.

### E. Provenance is substantially improved

The repository now contains runtime normalization, code hashes, mesh-quality records, replay checks, and inner-iteration histories. This strengthens reproducibility and makes it possible to distinguish a physical/numerical finding from a post-processing assumption.

## Main weaknesses and reviewer risks

### 1. RANS physical-step convergence is the principal blocker

This is the highest-risk issue. A reviewer can reasonably argue that the interpolation-method load difference is being compared below the demonstrated accuracy of the dual-time solve.

The clean solution is not another full 400-step matrix. A targeted control is enough if it shows stability.

### 2. Temporal convergence of the RANS response is not demonstrated

The RANS nondimensional step is 0.347224, much larger than the Euler value. The intended observable is a remeshing-induced transient; therefore time-step sensitivity matters even for a nominally steady physical problem.

A half-step control must replay the same physical mesh events rather than regenerate them from different iteration indices.

### 3. The RANS solver-control ladder is insufficiently documented

The files `rans_T04_sensitivity.csv` and `rans_T04_series.csv` are valuable, but labels L0/L1/L2 are not self-describing in the manuscript/repository summaries inspected here. Before they are used as evidence, archive the exact changed settings and the physical-event mapping for each level.

### 4. RANS spatial resolution is a limitation

The case has one viscous mesh class. This is acceptable for a paired transfer experiment if absolute aerodynamic validation is not claimed, but it prevents an absolute RANS accuracy statement. A reviewer at AIAA Journal may ask whether the remeshing signal is smaller than spatial-discretization uncertainty.

A full grid-convergence study for every method is unnecessary. One baseline mesh-sensitivity check, or a clear decision to avoid absolute-accuracy claims, is sufficient.

### 5. Random-sequence uncertainty is lightly sampled

The airfoil drift studies use two realizations; the Riemann stress matrix effectively uses one. This is enough to demonstrate large qualitative contrasts, but not to estimate a distribution or attach statistical confidence. The revised manuscript avoids such claims.

### 6. Native mesh mismatch remains inseparable in the target-referenced NACA norm

The second-/third-order target-referenced defect norms cluster tightly because the comparison includes the different native discrete equilibria. This means the norm should not be used to rank high-order transfer accuracy. The controlled history ratios are still valid and much stronger.

### 7. Boundary handling is now part of the method definition

GP2/GP3 projection followed by slip-wall enforcement is not the same mathematical operator as unconstrained L2 projection. The paper now states this, but the implementation description should consistently treat limiting, boundary projection, positivity floors, and vertex restriction as part of the returned transfer pipeline.

### 8. Cost is measured but not yet a strong engineering ranking

The RANS timing archive gives per-event transfer medians of roughly 5--7 s for GP2 and 3.5--5.5 s for NN2, while ten CFD steps take roughly 10--21 s depending on case and runtime conditions. These measurements show that transfer cost is not negligible, but system variability makes a clean method-cost ratio difficult. Keep timing secondary unless a production recommendation is made.

## Minimum work still required before submission

### P0-1 -- Close RANS physical-step convergence

Use one paired deformation sequence and a small representative subset:

- SA: NN2 and GP2;
- SST: NN2 and GP2;
- no-remeshing reference for each closure.

Either rerun a selected event window with a stricter/longer dual-time solve, or document and complete the existing L0/L1/L2 ladder if that already provides the required comparison.

**Acceptance criterion:** the reported remeshing-induced load quantities used in the paper (at least CD and CL, preferably pressure/friction split) must change by substantially less than the difference the manuscript attempts to resolve. If this criterion is not met, keep the RANS section but remove any comparative aerodynamic-accuracy conclusion.

### P0-2 -- Document the solver-control ladder

For L0/L1/L2 archive:

- exact TIME_STEP;
- INNER_ITER;
- convergence threshold and start iteration;
- linear-solver tolerances;
- number/timing of remeshing events;
- exact displacement files or hashes;
- initial histories.

Without this metadata the existing sensitivity CSV is diagnostic but not publication-grade verification.

### P0-3 -- Perform one RANS time-step check at fixed physical remeshing times

At minimum, halve the physical time step for one closure and both NN2/GP2 using the same physical deformation sequence. Double the step index between events so that the physical remeshing times are unchanged.

**Acceptance criterion:** the central qualitative conclusion (no resolved NN2/GP2 load ranking despite much smaller GP2 integral defects) must remain unchanged, and the magnitude of the reported response must be stable enough to quote.

### P0-4 -- Finalize the configuration/provenance archive

Before submission, store the exact configuration files and executable/source identifiers corresponding to every main case. The repository already records hashes and code provenance; make the exact runnable inputs available rather than relying only on tables in a report.

### P0-5 -- Rebuild the paper after the controls

Update only the numerical statements affected by P0-1/P0-3. Do not re-open the paper structure unless the controls reverse a conclusion.

## Recommended but not mandatory

### P1-1 -- One additional Riemann deformation realization

Repeat only the headline stress comparison (NN1, BAR2, GP2 at 400 transfers, SF4 or SF2) with one or two additional seeds. The current contrasts are so large that this is unlikely to change the qualitative conclusion, but it would remove an easy reviewer criticism.

### P1-2 -- One viscous-mesh sensitivity reference

Run the no-remeshing SA/SST baseline on one modestly refined mesh, or cite/compare a validated reference configuration using the same turbulence-model definitions. This is useful for AIAA Journal but not essential if the paper explicitly avoids absolute RANS validation claims.

### P1-3 -- Add cost only if making a practical recommendation

If the conclusions recommend GP2/NN2 operationally, normalize transfer wall time and extra nonlinear work on a controlled node allocation. Otherwise leave the current timing data as supplementary context.

## Work that is not required to close this paper

Do **not** add the following merely for completeness:

- pitching/plunging airfoils;
- moving-grid ALE/GCL analysis;
- a second geometry;
- three-dimensional cases;
- a full higher-order RANS matrix;
- another acoustic-radius sweep;
- many additional random seeds;
- direct Jacobian/adjoint solves.

These are natural follow-on studies but do not close the present evidentiary gap as efficiently as RANS numerical controls.

## Readiness by target journal

### AIAA Journal

**Not ready yet.** The paper has an aerospace-relevant contribution and a stronger scientific story than before, but the unresolved RANS physical-step/time-step sensitivity is likely to be a major reviewer concern because the paper discusses aerodynamic consequences. Complete the P0 controls first.

After those checks, AIAA Journal is a reasonable first target if the RANS load response is stable enough to support at least a bounded aerodynamic conclusion. If the controls show that method differences remain below numerical resolution, that negative result is still publishable if framed explicitly: the paper would then establish that transfer-integral improvements do not translate into a resolved load difference at the solver accuracy used.

### Computers & Fluids

**Close to ready after the same P0 controls.** The transfer-method, remapping, and temporal-history emphasis is an especially natural fit. The paper does not need a moving-body case for this venue.

### Journal of Computational Physics

**Not the preferred target in the current form.** The paper provides a careful comparative/diagnostic study rather than a sufficiently general new interpolation theory or algorithmic advance for the strongest JCP positioning.

## Final go/no-go criterion

I would call the paper **submission-ready** when all of the following are true:

- [ ] the RANS physical-step convergence check is documented and the retained load claims survive it;
- [ ] one time-step sensitivity comparison at fixed physical remeshing times is documented, or RANS aerodynamic claims are narrowed so that temporal resolution is no longer used to rank methods;
- [ ] L0/L1/L2 controls are fully defined or removed from the evidence chain;
- [ ] exact configurations/source hashes needed to reproduce the headline cases are archived;
- [ ] the final manuscript is rebuilt from the recalculated CSVs with no legacy result values remaining;
- [ ] figures/tables/captions are checked against the final CSVs and no result appears twice without a distinct purpose;
- [ ] the abstract and conclusions contain no stronger RANS statement than the numerical controls support.

The paper does **not** need more breadth. It needs the existing RANS evidence to be numerically closed.
