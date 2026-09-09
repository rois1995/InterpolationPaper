# Fixed NACA0012 BDF2 Restart Validation Protocol

## Objective

Use the converged fixed-airfoil calculations to validate the distinctive two-history BDF2 analysis of the paper, rather than only the spatial interpolation properties.

## Required native calculations

1. Converge an independent steady solution on source mesh A.
2. Converge the same physical problem independently on target mesh B.
3. Use identical physical models, reference quantities, force integration, and nonlinear tolerances.
4. Preserve both target states and source states in the exact restart format used by the solver.

The native target calculation is the reference for every transfer-induced quantity. Comparing only against the continuing source run mixes native mesh error with transfer error.

## Experiment 1: steady mixed-history decomposition

Let

- `UA_star` be the converged source state;
- `UB_star` be the converged target state;
- `Uhat = I_AB(UA_star)` be the source state transferred to B;
- `eps = Uhat - UB_star`.

Initialize four BDF2 restarts on target mesh B:

| case | history at n | history at n-1 | BDF2 forcing eta |
|---|---|---|---|
| native | `UB_star` | `UB_star` | `0` |
| full | `Uhat` | `Uhat` | `1.5 eps` |
| latest only | `Uhat` | `UB_star` | `2 eps` |
| older only | `UB_star` | `Uhat` | `-0.5 eps` |

Advance exactly one physical BDF2 step, fully converging the nonlinear/dual-time iterations. Use the same initial iterate for all four runs.

### Mandatory outputs

- componentwise mass-weighted norms of `eps` and `eta`;
- initial target BDF2 residual before nonlinear iterations;
- converged first-step state difference relative to the native target step;
- lift, drag, and pitching-moment differences;
- nonlinear iteration count and final residual;
- wall pressure coefficient and near-wall state differences.

### Direct checks

1. **Residual identity**

   For the same trial state, the restarted-minus-native BDF2 residual must equal

   `-(M_B eta)/dt`.

2. **Older-history test**

   The older-only case has zero latest-state error but must generally produce a nonzero first-step response. This directly demonstrates that the latest state alone is insufficient.

3. **Linear superposition and coefficient ratios**

   For sufficiently small transfer errors,

   `delta_full ≈ delta_latest + delta_older`,

   `delta_latest ≈ (4/3) delta_full`,

   `delta_older ≈ -(1/3) delta_full`.

   The same relations should approximately hold for each first-order aerodynamic-coefficient perturbation.

4. **Optional Jacobian test**

   If the target Jacobian/linear solver can be reused, solve

   `(1.5 M_B + dt J_B) delta_lin = M_B eta`

   and compare `delta_lin` with the actual converged first-step difference.

This experiment can be run immediately from converged steady solutions. It is the minimum direct validation of the paper's BDF2 history equation.

## Experiment 2: weakly unsteady fixed-airfoil extension

The steady experiment uses the same spatial defect in both history levels. To test the general case, retain the fixed airfoil but apply a small sinusoidal inflow-angle or vertical-gust perturbation:

`alpha_inf(t) = 8 deg + Delta_alpha sin(omega t)`.

A practical initial choice is `Delta_alpha = 0.25 deg` and reduced frequency `k = omega c/(2 U_inf) = 0.05`, provided the response remains approximately linear and is temporally resolved. Use at least 100 physical steps per forcing period, or demonstrate temporal convergence.

1. Run native A and native B under identical forcing until a repeatable periodic response is obtained.
2. Select at least eight phases over one period.
3. At each phase transfer both source histories to B.
4. Advance one step for the immediate-response analysis and 20-50 steps for short-time force propagation.
5. Repeat for all transfer methods; use raw and globally corrected NN/BAR variants at identical phases.

### Analysis

For each phase and method, compute

- `||eps_n||`, `||eps_nm1||`, and `||eta||`;
- `||delta_np1||`;
- initial history-residual norm;
- `Delta CL`, `Delta CD`, and `Delta Cm` at the first step and over the response window.

At each phase, compare the relationship of the first-step response with `||eta||` and with `||eps_n||`. The BDF2 model is supported when `eta` gives the tighter phase-conditioned relationship. Do not pool phases blindly because the target implicit operator and force sensitivity change with phase.

## Delayed-response extension

Create target meshes that are identical near the airfoil but differ in localized regions at increasing upstream distance. Transfer at the same physical state and record when the force defect first exceeds the native-target noise floor.

Report the support of the initial state defect and compare the measured onset with acoustic and convective travel-time estimates. This experiment should follow the immediate fixed-NACA validation; it is not required for the first steady batch.

## Recommended first batch

Run the steady mixed-history decomposition first for:

- NN1;
- NN2 or NN3;
- BAR2;
- BAR3;
- GP2;
- GP3;
- NN1+global correction;
- BAR2+global correction.

This gives a direct validation of the paper's central BDF2 relation with modest additional computational cost. The weakly unsteady phase sweep can then be restricted to the most informative methods if the full matrix is expensive.
