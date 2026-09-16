# Remaining simulation plan for full closure of the InterpolationPaper

**Repository basis:** `rois1995/InterpolationPaper`, release `v1.0` / commit `6726ba2968fd8a32fc4e6517588c18be91db390f`

## Purpose

The current manuscript is already defensible for submission because its claims have been narrowed to what the existing controls support. The simulations below are therefore **not all required for first submission**. They are organized according to how much scientific uncertainty they would remove.

The key distinction is:

- **Tier A — close explicit remaining gaps in the current paper.** These are compact, high-value runs.
- **Tier B — close methodological caveats that are currently scoped honestly.** Useful if time is available.
- **Tier C — broaden the claims.** Not recommended before submission unless the paper is deliberately expanded.

No moving-body, pitching/plunging, three-dimensional, or additional-geometry campaign is required to close the present fixed-domain paper.

---

# Tier A — High-value closure runs

## A1. SST time-step replay at half the production physical step

### Gap
The current temporal-sensitivity test is available only for SA. At the production RANS time step, the paired NN2–GP2 aerodynamic difference is smaller than the time-step sensitivity measured in SA. The manuscript therefore correctly refuses to rank NN2 and GP2 aerodynamically, but SST has no equivalent direct temporal-control calculation.

### Runs
Use the existing SST production case and the **same three physical remeshing events** already used in the control window.

Run:

1. `SST / no-transfer reference / Δt/2`
2. `SST / NN2 / Δt/2`
3. `SST / GP2 / Δt/2`

The production controls span physical steps 2–40 with events corresponding to production event indices 9, 19, and 29. At half step, replay the exact archived displacement files at the doubled time indices, exactly as already done for the SA half-step control.

Do **not** regenerate the random displacements using the new step index. Reuse the archived displacement arrays by physical event.

### Keep fixed
- geometry and mesh;
- turbulence model and freestream;
- transfer settings;
- three physical event times;
- displacement arrays;
- force/reference definitions;
- convergence criterion;
- all numerical settings except physical `TIME_STEP`.

### Compare
For `CL`, `CD`, `CDp`, `CDv`:

- case-minus-reference response for NN2;
- case-minus-reference response for GP2;
- paired difference
  \[
  P(t)=D_{\mathrm{NN2}}(t)-D_{\mathrm{GP2}}(t);
  \]
- RMS and maximum change relative to production;
- endpoint change only as a secondary diagnostic.

### Interpretation
If the SST paired difference changes by an amount comparable to or larger than its production RMS, the existing conclusion is confirmed: the aerodynamic method difference is not temporally resolved.

If the paired difference is stable, SST can be discussed more strongly than SA, but a **third time-step level** would still be required before claiming formal temporal convergence.

### Manuscript consequence
This closes the present asymmetry between SA and SST temporal controls.

### Cost
Low: three short control trajectories, not full 400-step production reruns.

---

## A2. SST inner-convergence closure beyond the 1600-iteration cap

### Gap
The cap-1600 SST controls still terminate at the cap for 28 of 38 analyzed physical steps. Their load changes are small, which is reassuring, but this is a sensitivity test rather than a fully converged dual-time reference.

### Runs
For the same SST seed-1 control window and same three remeshing events:

1. `SST / NN2 / higher-inner-convergence control`
2. `SST / GP2 / higher-inner-convergence control`

The no-transfer SST reference is already inexpensive and well converged, so a new reference run is needed only if its numerical settings must be changed for strict consistency.

Use one of:

- `INNER_ITER = 3200`, retaining the same residual target; or preferably
- a configuration that continues until the `log10(rms rho) = -12` criterion is actually met, with a sufficiently high safety cap.

### Required logging
For every physical step archive:

- inner iterations used;
- final density residual;
- turbulence residuals;
- force coefficients at physical-step convergence;
- whether the cap was reached.

### Acceptance criterion
For `CL` and `CD`, require the change from the 1600-control result to be small compared with the production NN2–GP2 paired difference. A practical target is:

- individual case-minus-reference RMS change < 1–2% of the response RMS;
- paired NN2–GP2 RMS change < 10% of the paired difference.

If the pair remains more sensitive than this, retain the current conservative wording and do not attempt a RANS aerodynamic ranking.

### Manuscript consequence
Would allow the paper to replace “cap sensitivity is small” with a stronger statement that iterative error is demonstrably below the compared aerodynamic signal.

### Cost
Very low to low: two short 38-step controls.

---

## A3. Replicate the NN1 + global-correction Riemann stress test for seeds 2 and 3

### Gap
The raw NN1, BAR2, and GP2 broadening hierarchy is now reproduced across three deformation realizations. The strongest conceptual statement—global correction restores the density integral but does not recover the spatial profile—is still demonstrated for only the seed-1 sequence.

### Runs
At 400 transfers, run NN1 with global correction for:

1. SF4, seed 2
2. SF4, seed 3
3. SF2, seed 2
4. SF2, seed 3

Use exactly the same solver settings and profile extraction currently used for `riemann_amplitude_seeds.csv`.

### Metrics
At step 399:

- position error in cells;
- 10–90% width change in cells;
- cumulative interpolation-only density integral change;
- density/energy integral defect per transfer if already available.

### Expected diagnostic
The key test is not that every realization gives the same width. It is whether:

\[
|\Delta I_\rho| \rightarrow O(\epsilon_{\mathrm{mach}})
\]

while the profile remains substantially broader than the no-remeshing reference.

### Acceptance criterion
For each new realization:

- corrected cumulative density integral near roundoff;
- width change remains materially larger than GP2 and materially nonzero.

### Manuscript consequence
This would turn the global-conservation counterexample from a single-realization result into a replicated result. Given the importance of this claim and the low cost of the Euler test, this is probably the **highest value/cost simulation still available**.

### Cost
Low: four inexpensive Euler/Riemann runs.

---

# Tier B — Close currently scoped methodological caveats

## B1. BDF2 zero-forcing cancellation experiment

### Gap
The full/latest/older experiment validates the BDF2 coefficients along one spatial discrepancy direction. It does not directly demonstrate that two individually nonzero history perturbations cancel when their BDF2 combination vanishes.

### Construction
On the existing Euler NACA target mesh, use a boundary-admissible discrepancy vector \(\epsilon\) and define

\[
\epsilon^n=\frac14\epsilon,\qquad
\epsilon^{n-1}=\epsilon.
\]

Then

\[
\eta^n
=
2\epsilon^n-\frac12\epsilon^{n-1}
=
0.
\]

### Runs
For one pointwise method (BAR2 or NN2) and optionally GP2 after explicitly enforcing the same slip-wall admissibility used by the receiving solver:

- native history pair;
- cancellation history pair.

Only the first converged BDF2 physical step is required.

### Expected result
The two solved states should agree to the established nonlinear/iterative floor because the history forcing is zero.

### Metrics
- volume-weighted state norm;
- `ΔCL`, `ΔCD`, `ΔCMz`;
- common-trial residual difference before solve;
- wall/interior split if GP2 is included.

### Manuscript consequence
Closes the most direct algebraic caveat left by the collinear full/latest/older experiment.

### Cost
Negligible: one-step restart calculations.

---

## B2. Perturbation-amplitude sequence for the BDF2 linear-response approximation

### Gap
The measured state superposition defect is approximately 1.8–3.5%. The current paper correctly calls this local directional evidence rather than validation of the full linearized resolvent.

### Runs
Choose one representative pointwise method and GP2, preferably BAR2 + GP2.

Scale the same target-referenced discrepancy by:

\[
\lambda=1,\quad \frac12,\quad \frac14.
\]

For each \(\lambda\), run:

- full history;
- latest-only;
- older-only.

This is 9 one-step runs per method if all configurations are repeated; the native target step is reusable.

### Metrics
For each amplitude:

\[
e_{\mathrm{sup}}^U
=
\frac{
\|\delta_{\mathrm{full}}
-\delta_n-\delta_{n-1}\|
}{
\|\delta_{\mathrm{full}}\|
},
\]

plus the corresponding lift superposition error.

### Expected result
For a smooth local nonlinear response, the absolute nonlinear remainder should decrease faster than the first-order response. Consequently the **normalized** superposition defect should decrease approximately linearly with perturbation amplitude until it reaches the solver floor.

### Manuscript consequence
Would upgrade “directional consistency with the BDF2 coefficients” to a genuine perturbation-amplitude verification of the local linearization.

### Cost
Very low.

---

## B3. Fine-mesh remeshing control window

### Gap
The current finer RANS mesh is used only for a no-transfer load sensitivity. It shows that absolute load values remain spatially sensitive, but it does not test whether the *remeshing response* or the NN2–GP2 difference is mesh-sensitive.

### Runs
On the existing finer viscous mesh, for one closure first (SA is sufficient as the primary test):

1. fine-mesh no-transfer reference;
2. fine-mesh NN2, same three physical remeshing events;
3. fine-mesh GP2, same three physical remeshing events.

If the result is scientifically important, repeat with SST.

### Important design requirement
The random perturbation must represent the same **physical remeshing severity**, not merely the same node-wise random numbers. Match:

- physical displacement amplitude relative to local mesh size;
- protected wall region;
- event times;
- quality threshold.

### Compare
- case-minus-reference `CL`, `CD`, `CDp`, `CDv`;
- paired NN2–GP2 difference;
- transfer integral defects;
- inner-iteration effort.

### Manuscript consequence
Would determine whether the application-level method comparison is stable to mesh refinement. This is more relevant to the paper than a large formal drag-validation exercise.

### Cost
Moderate but still much smaller than repeating the complete 400-step matrix.

---

# Tier C — Only if broader claims are desired

## C1. Formal three-level RANS spatial-convergence study

The current coarse/fine comparison is a **sensitivity check**, not a formal uncertainty estimate.

A true spatial-convergence study requires a systematically related third mesh and controlled refinement ratios. Ideally refine:

- surface spacing;
- wake resolution;
- outer triangles;
- wall-normal layering/first-cell spacing in a controlled way.

Run steady/no-transfer SA and SST on all three levels.

This is needed only if the paper is to report formal spatial convergence, GCI-like uncertainty, or absolute RANS aerodynamic accuracy. It is **not required** for the current transfer/remeshing claims.

---

## C2. Third physical-time level

If formal temporal convergence of the RANS remeshing transient is desired, use

\[
\Delta t,\quad \Delta t/2,\quad \Delta t/4
\]

with identical physical remeshing times and replayed displacement files.

This should be done for a **short control window**, not for the full 400-step campaign.

Run reference + NN2 + GP2 for the selected closure. SA is the natural first target; repeat SST only if a closure-independent temporal-order claim is desired.

This is required only for a formal temporal-order/uncertainty statement.

---

## C3. Additional random realizations for the RANS airfoil

The present RANS study has two deformation realizations. That is sufficient for the paper's current non-statistical conclusions but insufficient to characterize a seed distribution.

If you want statistically stronger statements, add at least three more common sequences:

- SA NN2/GP2;
- SST NN2/GP2.

This means 12 additional full trajectories for three new seeds.

Do this **only after** temporal sensitivity is under control; otherwise more seeds quantify variability of a numerically under-resolved observable.

---

# Analysis-only closure items — no new CFD simulation required

## D1. SST post-floor integral accounting
The archived SST transfer integral is evaluated before the primitive-variable positivity floor. Using the stored pre-floor and returned restart states, compute the integral:

1. before transfer;
2. after raw transfer;
3. after positivity floor;
4. after any boundary/admissibility operation.

This would define the conservation properties of the **actual returned restart state**.

## D2. Transfer timing/cost
The timing data already exist. If an engineering recommendation is added, normalize:

- transfer wall time per event;
- mesh-motion time;
- CFD time between events;
- additional inner iterations caused by remeshing.

No new run is required unless timing reproducibility is considered important.

## D3. Boundary-admissibility audit for GP
The Euler wall-normal momentum mechanism is already identified. A compact post-processing table showing raw GP transfer → boundary-admissible history → solver residual would make the result easier to audit without additional simulations.

---

# Recommended execution order

| Priority | Task | New CFD runs | Scientific value | Required for current submission? |
|---|---|---:|---|---|
| **1** | NN1+GC Riemann seeds 2/3 at SF4/SF2 | 4 cheap Euler runs | Very high | No, but strongly recommended |
| **2** | SST half-step replay, reference + NN2 + GP2 | 3 short RANS runs | Very high | No |
| **3** | SST higher-inner-convergence NN2/GP2 | 2 short RANS runs | High | No |
| **4** | BDF2 cancellation pair | 1–2 one-step runs | High / very cheap | No |
| **5** | BDF2 amplitude sequence | ~18 one-step runs for 2 methods | High / cheap | No |
| **6** | Fine-mesh remeshing window | 3 SA runs; +3 SST if desired | Medium-high | No |
| **7** | Third spatial mesh | several steady/reference runs | Medium | Only for formal spatial uncertainty |
| **8** | Δt/4 temporal level | 3 short runs per closure | Medium | Only for formal temporal convergence |
| **9** | Additional RANS seeds | 12 full runs for +3 seeds | Low before convergence closure | No |

---

# Recommended stopping point

If the aim is to **submit the current paper without expanding its scope**, I would stop after A1–A3, and optionally B1 because it is exceptionally cheap.

That gives:

- replicated evidence that global conservation does not repair spatial diffusion;
- symmetric SA/SST time-step controls;
- a genuinely converged SST inner-solve check;
- optionally an exact BDF2 zero-forcing cancellation test.

After those runs, there would be very little left for a reviewer to attack without asking for a fundamentally broader paper.

If the aim is instead to make a formal quantitative statement about **which interpolation method gives more accurate RANS aerodynamic loads**, then B3 + C1 + C2 become necessary. That is a materially larger project and is not needed for the paper's present contribution.
