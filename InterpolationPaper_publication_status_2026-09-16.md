# Publication-readiness assessment — 2026-09-16

## Executive decision

**The manuscript is ready for journal submission after the accompanying wording and presentation pass.** I do not consider any additional simulation campaign mandatory before first submission, provided the RANS claims remain scoped as they are in the revised manuscript: the RANS calculations demonstrate transfer-state/integral behavior and bound the aerodynamic comparison, but they do not establish a quantitative NN2-versus-GP2 load-accuracy ranking.

This decision is stronger than the previous readiness assessment because the repository now contains the controls that were previously missing: three independent Riemann deformation realizations, executed-run runtime normalization, an increased dual-time iteration-cap control for SA and SST, a half-step replay control for SA at the same physical remeshing times, a finer-mesh no-transfer sensitivity check, code/checkpoint provenance, and a clean release/regeneration workflow.

## Changes made in this manuscript pass

1. **RANS inner-convergence wording narrowed.** The 1600-iteration control changes the paired lift/drag difference by roughly 6--14%, but 28 of 38 SST control steps still reach the larger cap. The manuscript therefore no longer calls the SST load response fully inner-converged. Instead, it states that the 800-iteration cap is not the dominant observed sensitivity for lift/drag while the SST control remains an incomplete convergence proof.
2. **Temporal-control scope made explicit.** The half-step replay exists only for SA. The revised text no longer generalizes its quantitative result to SST.
3. **Mesh refinement interpreted as sensitivity, not uncertainty.** The finer viscous mesh provides one additional resolution level. It demonstrates appreciable absolute-load mesh sensitivity but cannot support a formal discretization-uncertainty estimate or grid-convergence claim.
4. **RANS conclusion sharpened.** GP2 clearly returns smaller monitored integral defects, but the aerodynamic NN2--GP2 difference is not time-step resolved in the SA control. The paper therefore reports the absence of a resolved load ranking rather than aerodynamic equivalence.
5. **Riemann replication added to the main argument.** The broadening hierarchy is reproduced over three deformation realizations at both deformation amplitudes, materially strengthening the repeated-remeshing conclusion.
6. **Appendix presentation improved.** The numerical-control table is restricted to lift and total drag, float placement is controlled, and captions state exactly what each control does and does not establish.
7. **Data-provenance language made publication-facing.** Internal development/audit language was removed; the text now describes the versioned source tables, generated manuscript views, replay controls, and release checks directly.

## Evidence strength

### 1. Smooth interpolation verification — strong

The implemented NN, barycentric, and Galerkin variants recover the expected smooth-field orders on the tested structured and mixed mesh families. This is direct code-verification evidence for the transfer operators and is independent of the CFD solver response.

### 2. Boundedness versus conservation — strong

The discontinuous test cleanly separates nominal order, boundedness, and conservation. Unlimited higher-order reconstructions can create extrema; limiting can restore bounds while worsening the lumped integral; global correction can restore the monitored integral without determining the spatial error. This is a strong supporting block because each property is measured independently.

### 3. Repeated Riemann remeshing — strong and now reproducible

The common tracker identifies cumulative broadening as the dominant low-order repeated-transfer error. The large NN1/BAR2 broadening and very small GP2 broadening are reproduced across three deformation realizations. Global correction removes the NN1 integral defect to roundoff without removing the broadening. This is one of the paper's cleanest counterexamples to using global conservation as a spatial-accuracy criterion.

### 4. Controlled BDF2 history experiment — strongest distinctive result

The older-only case gives a nonzero state and load response with an exact latest history, directly demonstrating that a latest-state-only restart diagnostic is insufficient for BDF2. The directional response ratios recover the BDF2 coefficients closely, and the residual audit explains the Galerkin wall-normal discrepancy through slip-wall admissibility. This is the clearest contribution beyond a conventional interpolation-method comparison.

### 5. Euler aerodynamic response — useful, configuration-specific

The repeated Euler NACA case gives a reproducible GP2/NN2 drag-trend separation over two deformation realizations, while the localized-transfer experiment supports an acoustic-scale delayed response. These results are useful but should remain explicitly configuration-specific; they are not evidence of a universal method ranking.

### 6. Fixed-body RANS extension — sufficient as a scoped application

The RANS extension now has enough numerical controls to be publishable **as a bounded application result**:

- GP2 reduces monitored turbulence-inventory defects by factors of about 37--247 relative to NN2.
- Increasing the dual-time cap changes the paired lift/drag difference by only about 6--14%, although SST still reaches the larger cap frequently.
- In SA, halving the physical step while replaying the same remeshing events changes the time-aligned paired method difference by more than its production RMS. Hence the paired transient is not temporally converged.
- A single finer mesh changes the absolute steady loads materially; this is mesh sensitivity, not a formal uncertainty estimate.

These controls are sufficient to justify the paper's actual conclusion: **a much smaller transfer-integral defect does not, at the present numerical resolution, produce a resolved aerodynamic advantage.** They are not sufficient to rank NN2 and GP2 by RANS load accuracy, and the manuscript does not do so.

## Remaining weaknesses and likely reviewer questions

### RANS temporal convergence

This remains the principal technical limitation, but it is now measured rather than unknown. The SA half-step control shows that the detailed paired transient is not time-step converged. There is no SST half-step control. A reviewer may request one, but it is no longer required to support the manuscript's deliberately limited RANS conclusion.

### SST inner convergence

The 1600-iteration SST control still reaches its cap in 28 of 38 tested steps. The load changes are modest relative to the transfer response and method-difference scale, but this is not a complete asymptotic convergence demonstration. The manuscript states this directly.

### Spatial convergence of RANS loads

Only one finer mesh is available, with the same first off-wall spacing. It establishes sensitivity, not an order or uncertainty. This is acceptable because the paper does not claim RANS validation or absolute aerodynamic accuracy.

### Statistical breadth

The Euler airfoil comparison uses two deformation realizations; the replicated Riemann stress test uses three. These are sufficient for the qualitative and paired conclusions reported, but not for population-level statistical claims. None are made.

### Scope

The study is two-dimensional and uses fixed physical domains with discrete mesh replacement. Moving meshes, pitching/plunging bodies, dynamic stall, and 3-D applications are outside the stated scope. They should remain future work rather than submission requirements.

### Novelty positioning

The interpolation operators themselves are established. The paper's contribution is the controlled decomposition of remeshing error in a multistep CFD pipeline, particularly the two-history BDF2 experiment, the conservation-versus-spatial-error counterexamples, boundary-admissibility audit, and the demonstration that transfer diagnostics do not map monotonically to aerodynamic outputs. The introduction and conclusions should continue to position novelty there rather than in the existence of NN, barycentric, or Galerkin transfer.

## Submission status by venue

### AIAA Journal

**Ready for submission with the revised claims.** The aerospace relevance is supplied by the aerodynamic restart/load experiments and the RANS application. The most likely reviewer challenge is the RANS temporal-control limitation. The paper already quantifies that limitation and does not base its main contribution on a resolved RANS method ranking, so I would not delay first submission solely to add another campaign.

### Computers & Fluids

**Ready for submission.** The balance of transfer-method verification, conservation/boundedness analysis, repeated-remeshing behavior, and multistep restart diagnostics fits a numerical-CFD audience particularly naturally.

## Optional strengthening before or during review

These would improve the paper but are not prerequisites for first submission:

1. Run the SST half-step replay on the same physical remeshing events, mainly to parallel the SA temporal control.
2. Add a third viscous mesh only if a reviewer requests a formal RANS grid-convergence statement; otherwise retain the current sensitivity wording.
3. Archive the tagged release in a persistent repository such as Zenodo and cite the DOI in the data-availability statement.
4. Add cost comparisons only if the manuscript is expanded toward a production-method recommendation.

## Go / no-go

**GO for submission after this revision.**

The central scientific claims are supported by direct, reproducible evidence. The remaining limitations affect claims the manuscript explicitly declines to make, rather than invalidating its stated contribution.
