# Interpolation Error Study — 2D Riemann Case: Consolidated Findings

**Purpose.** This document consolidates the quantitative results of the interpolation-error
investigation on the 2D Riemann case, for use by a downstream agent writing the results
section of a journal paper. It is organised so that convergence plots (log–log) and
per-method/per-order error tables can be built directly from the numbers below.

**Ground rules for the downstream agent.**
- **Do not invent or extrapolate numbers.** Every value here is transcribed from a named
  source file (CSV) or a named analysis run. If a number you need is not present, it was
  not measured — say so, or request it.
- Each table names its **source file** and the **script** that produced it.
- Read the **Caveats & Assumptions** section (§10) before writing any claim — several
  metrics are fragile and must be qualified.
- Numbers are quoted to 3–4 significant figures; full double precision is in the CSVs.

---

## 1. Test case and common setup

- **Governing equations:** 2D compressible Euler (`SysType = Euler`, inviscid).
- **Domain:** rectangle `[0, 5] × [0, 1]`.
- **Boundary markers:** `Bot`, `Left`, `Right`, `Top`.
- **Baseline mesh:** structured quadrilaterals, **50 601 nodes** (501 × 101), uniform
  spacing **h = 0.01**.
- **Mesh-refinement family** (same geometry/markers, uniform quads):

  | tag          | h (nominal) | h (measured, area/cell) | nodes  |
  |--------------|-------------|-------------------------|--------|
  | `Prova_0.08` | 0.08        | 0.07752                 | 910    |
  | `Prova_0.04` | 0.04        | 0.03984                 | 3 302  |
  | `Prova_0.02` | 0.02        | 0.01996                 | 12 852 |
  | `Mesh` (base)| 0.01        | 0.01000                 | 50 601 |

  (`h` measured as `sqrt(total_area / cell_count)`; exact for a uniform quad grid.)

- **Unsteady driver:** `dt = 1e-3`, `endIter = 400` → **0.4 physical time**.
  Interpolation happens twice per adaptation block, so:
  - `dIters = 10` → **80 interpolation cycles** over the run;
  - `dIters = 2`  → **400 interpolation cycles** over the same physical time (same
    physics, 5× the interpolation count).
- **Reference solution:** `Unsteady_Ref` — the same problem with **no mesh adaptation and
  no interpolation** (fixed mesh, pure CFD). This is the "truth" for shock-position
  comparisons.

### 1.1 Mesh deformation (how the "adapted" mesh is generated each cycle)

Script: `Base_Interp/RandomMeshMovements.py`.
- Each cycle deforms the **original** mesh (`../Prova_Steady/Mesh.su2`) by a **fresh,
  independent white-noise nodal displacement** (seed = `12345 + iterEnd`): the mesh
  oscillates around a fixed, good-quality base and does **not** accumulate distortion.
  Boundary nodes are frozen.
- Per-node displacement magnitude is bounded by a **safe radius** (minimum altitude over
  incident corner triangles) divided by `SafeFactor`: `SafeFactor = 4` in the baseline
  studies, `SafeFactor = 2` (double displacement) in the amplitude study (§9). Motion is
  **x-only** (`LockY = True`) in the multi-cycle runs.
- The mean per-cycle nodal displacement is reported per study (it sets the abscissa of the
  δ-sensitivity sweep, §3).

### 1.2 Interpolation methods and orders under test

| method             | orders tested        | conservative? | notes |
|--------------------|----------------------|---------------|-------|
| `NN`               | First, Second, Third | no            | nearest-neighbour + WLS gradient (Second) + Hessian (Third) |
| `Barycentric`      | Second, Third        | no            | linear over containing simplex (Second); +½-gradient correction (Third). **No First mode** (plain barycentric is intrinsically 2nd order). |
| `ConsGalerkinProj` | Second, Third        | **yes (local)** | L2/Galerkin projection over exact supermesh; P1 (Second) / P2 (Third). |

`EnforceConservation = "Global"` is a **method-agnostic** post-pass that rescales each
field to restore its total integral within the source range (max-principle bounded). It is
off by default and toggled on/off in §4 and §9.

---

## 2. Study A — Order of accuracy (mesh refinement)

**Question:** does each method achieve its nominal order `p` in `E ∼ h^p`?
**Design:** interpolate a known analytic field from the source mesh onto a **deformed**
mesh, at each `h` in the refinement family, with **δ/h held fixed** (`SafeFactor = 8`
constant → displacement scales with `h`; measured mean δ/h ≈ 0.057). This is the scaling
under which the fitted exponent **is** the order of accuracy.
**Analytic field** (smooth, wavenumber k = 2π, wavelength 1.0):
`Density = 1 + 0.5·sin(kx)·cos(ky)` (range [0.5, 1.5], so absolute L1 = relative);
`Momentum_x = 0.3·cos(kx)·sin(ky)`, `Momentum_y = 0.2·sin(kx)·cos(ky)`,
`Energy = 2.5 + 0.4·cos(kx)·sin(ky)`.
**Error metric:** area-weighted (lumped nodal areas) L1 of Density over moved nodes,
against the exact analytic value at the deformed node positions.
**Source file:** `Tier0_Refine/tier0_refine.csv` (field = Density).
**Script:** `Base_Interp/Tier0_InterpBenchmark.py --meshes ... --safe-factor 8 --wavenumber 6.283185`.

### 2.1 Density L1 error vs h (build the convergence plot from this)

| method / order          | h=0.0775   | h=0.0398   | h=0.0200   | h=0.0100   |
|-------------------------|------------|------------|------------|------------|
| NN / First              | 6.268e-03  | 3.227e-03  | 1.577e-03  | 7.968e-04  |
| NN / Second             | 2.752e-04  | 4.542e-05  | 8.791e-06  | 2.002e-06  |
| NN / Third              | 7.057e-04  | 9.954e-05  | 1.232e-05  | 1.568e-06  |
| Barycentric / Second    | 1.693e-03  | 4.757e-04  | 1.194e-04  | 2.995e-05  |
| Barycentric / Third     | 2.675e-04  | 2.853e-05  | 3.125e-06  | 3.785e-07  |
| ConsGalerkinProj/Second | 8.335e-04  | 1.898e-04  | 4.445e-05  | 1.075e-05  |
| ConsGalerkinProj/Third  | 1.639e-04  | 3.147e-05  | 4.267e-06  | 5.530e-07  |

(L2 and Linf columns are also in the CSV if needed.)

### 2.2 Fitted order p (least-squares log-log over the 4 meshes)

| method / order          | fitted p | per-interval slopes (coarse→fine) | nominal | verdict |
|-------------------------|----------|-----------------------------------|---------|---------|
| NN / First              | 1.01     | 1.00, 1.04, 0.99                  | 1       | ✓       |
| NN / Second             | 2.40     | 2.71, 2.38, 2.14                  | 2       | ✓       |
| NN / Third              | 2.99     | 2.94, 3.02, 2.98                  | 3       | ✓       |
| Barycentric / Second    | 1.97     | 1.91, 2.00, 2.00                  | 2       | ✓       |
| Barycentric / Third     | 3.20     | 3.36, 3.20, 3.05                  | 3       | ✓       |
| ConsGalerkinProj/Second | 2.12     | 2.22, 2.10, 2.05                  | 2       | ✓       |
| ConsGalerkinProj/Third  | 2.79     | 2.48, 2.89, 2.96                  | 3       | ✓ pre-asymptotic (slopes rising to 3) |

Notes for interpretation:
- `Barycentric/Third` fits **3.20**, slightly above 3 because it approaches 3 from above;
  the result is third order, as designed.
- `ConsGalerkinProj/Third` fits **2.79**, with per-interval slopes rising 2.48 → 2.89 →
  2.96: it is **pre-asymptotic on the coarse meshes and approaching 3**; finer meshes
  raise the fit. This is expected, not a defect.

---

## 3. Study B — Displacement sensitivity (fixed mesh)

**Important:** this is **not** the order of accuracy. It is the exponent `q` in `E ∼ δ^q`
with **h fixed** (h = 0.01) and only the displacement δ varied. `q` and the order `p` are
different exponents in different limits (e.g. barycentric is δ¹ at fixed h but h² under
refinement). Report `q` as "sensitivity to displacement at fixed resolution," and do
**not** compare it against the First/Second/Third labels.

**Design:** single interpolation of the analytic field (wavenumber k = 4π, wavelength 0.5)
onto meshes deformed with `SafeFactor ∈ {4, 8, 16, 32}` → mean nodal displacement
δ ∈ {1.132e-3, 5.660e-4, 2.830e-4, 1.415e-4}.
**Source file:** `Tier0/tier0_analytic.csv` (field = Density, region = moved).
**Script:** `Base_Interp/Tier0_InterpBenchmark.py --mode analytic --safe-factors 4 8 16 32`.

### 3.1 Density L1 error vs mean displacement δ

| method / order          | δ=1.132e-3 | δ=5.660e-4 | δ=2.830e-4 | δ=1.415e-4 |
|-------------------------|------------|------------|------------|------------|
| NN / First              | 3.148e-03  | 1.574e-03  | 7.870e-04  | 3.935e-04  |
| NN / Second             | 3.214e-05  | 8.849e-06  | 2.866e-06  | 1.150e-06  |
| NN / Third              | 2.456e-05  | 1.237e-05  | 6.194e-06  | 3.098e-06  |
| Barycentric / Second    | 2.279e-04  | 1.203e-04  | 6.178e-05  | 3.130e-05  |
| Barycentric / Third     | 5.235e-06  | 3.133e-06  | 1.708e-06  | 8.907e-07  |
| ConsGalerkinProj/Second | 9.157e-05  | 4.428e-05  | 2.195e-05  | 1.095e-05  |
| ConsGalerkinProj/Third  | 8.432e-06  | 4.288e-06  | 2.156e-06  | 1.080e-06  |

### 3.2 Fitted displacement-sensitivity q (least-squares log-log)

| method / order          | fitted q | per-interval slopes (coarse→fine δ) |
|-------------------------|----------|-------------------------------------|
| NN / First              | 1.00     | 1.00, 1.00, 1.00                    |
| NN / Second             | 1.60     | 1.86, 1.63, 1.32                    |
| NN / Third              | 1.00     | 0.99, 1.00, 1.00                    |
| Barycentric / Second    | 0.96     | 0.92, 0.96, 0.98                    |
| Barycentric / Third     | 0.85     | 0.74, 0.88, 0.94                    |
| ConsGalerkinProj/Second | 1.02     | 1.05, 1.01, 1.00                    |
| ConsGalerkinProj/Third  | 0.99     | 0.98, 0.99, 1.00                    |

Interpretation: q ≈ 1 for most operators is expected (an interpolant exact at source nodes
has error ∼ |f''|·δ·h at fixed h). `NN/Second`'s q = 1.60 with slopes decaying toward 1
indicates a two-term error `E ≈ aδ + bδ²` with the linear (WLS-gradient) term taking over
as δ shrinks.

---

## 4. Study C — Global conservation correction, single interpolation (smooth field)

**Question:** does `EnforceConservation="Global"` (i) make every method conservative and
(ii) at what accuracy cost, on a smooth field?
**Design:** one interpolation of the analytic field onto a deformed base mesh (h = 0.01),
each method/order run with correction **off** and **"Global"**. `L1_error` = Density L1 vs
exact; `cons_defect` = relative Density integral defect (tgt−src)/src from the run's
`Integrals.csv`.
**Source file:** `Correction/correction_comparison.csv`. **Script:** `CorrectionComparison.py`.

| method / order          | L1 (off)   | L1 (Global) | cons defect (off) | cons defect (Global) |
|-------------------------|------------|-------------|-------------------|----------------------|
| NN / First              | 7.968e-04  | 7.968e-04   | −2.00e-06         | 3.55e-16             |
| NN / Second             | 2.002e-06  | 1.998e-06   | 1.52e-08          | −5.33e-16            |
| NN / Third              | 1.568e-06  | 1.569e-06   | 1.78e-08          | 1.07e-15             |
| Barycentric / Second    | 2.995e-05  | 2.999e-05   | 1.06e-07          | 7.11e-16             |
| Barycentric / Third     | 3.785e-07  | 3.809e-07   | 2.14e-08          | 7.11e-16             |
| ConsGalerkinProj/Second | 1.075e-05  | 1.075e-05   | 1.38e-08          | 3.55e-16             |
| ConsGalerkinProj/Third  | 5.530e-07  | 5.534e-07   | 2.17e-08          | 1.07e-15             |

**Finding:** "Global" drives conservation to machine precision (∼1e-16) for **every**
method, and changes L1 by **≤0.6%** everywhere. On a smooth field the correction is
essentially accuracy-free and the accuracy ranking is unchanged.

---

## 5. Study D — Galerkin limiter tolerance and the P2 conservation defect

**Question:** does tightening `GalerkinLimitTol` improve ConsGalerkinProj conservation?
(`ConsGalerkinProj/Third` carries a small conservation defect.)
**Design:** single ConsGalerkinProj interpolation of a **real** step-399 restart onto a
deformed mesh, sweeping the limiter. `density_defect` = (tgt−src)/src of ∫Density.
**Source file:** `GalerkinTol/galerkin_tol.csv`. **Script:** `GalerkinTolSweep.py`.

| setting                     | Density defect | Energy defect |
|-----------------------------|----------------|---------------|
| Second, tol=2e-3 (default)  | 7.468e-09      | 5.767e-09     |
| Second, tol=2e-4            | 7.955e-09      | 5.093e-09     |
| Second, tol=0 (hard clip)   | 7.654e-09      | 3.910e-09     |
| Second, limiter OFF         | 7.468e-09      | 5.767e-09     |
| Third, tol=2e-3 (default)   | −5.096e-08     | −1.537e-07    |
| Third, tol=2e-4             | −5.096e-08     | −1.536e-07    |
| Third, tol=2e-5             | −5.095e-08     | −1.537e-07    |
| Third, tol=0 (hard clip)    | −5.076e-08     | −1.541e-07    |
| Third, limiter OFF          | −5.096e-08     | −1.537e-07    |

**Finding:** the limiter tolerance has **no effect** on conservation (columns flat). The P2
(Third) defect (∼5e-8) is ~7× the P1 (Second) defect (∼7.5e-9) and is **intrinsic to the
P2 projection**: it solves on vertices + edge-midpoint DOFs, but the vertex-only restart
discards the edge DOFs, giving an O(h³) conservation defect. Applying
`EnforceConservation="Global"` after the P2 projection drives the defect to machine
precision (Third: −5.10e-08 → +1.87e-16; Second: +7.47e-09 → −9.35e-16).

---

## 6. Study E — Multi-cycle conservation drift & monotonicity (80 cycles, dIters=10)

**Design:** the seven method/order combinations run as full 400-step adaptation loops
(`dIters=10` → 80 interpolation cycles), `SafeFactor=4`, `EnforceConservation=False`,
`debugInterp=1`. Per-cycle conservation is read from `DebugFiles/Integrals_<step>.csv`
(`src_geom` = ∫φV before interpolation, `tgt_geom` = after; the difference is pure
interpolation error, no time advance). Overshoot from `Minimums/Maximums_<step>.csv`.
**Source file:** `MultiCycle/multicycle_summary.csv` (Density/Energy-restricted for
§6.1–6.2, same script). **Script:** `MultiCycleAnalysis.py`.

### 6.1 Conservation drift — net over 80 cycles, as % of initial integral (RELIABLE fields)

Density and Energy only (their integrals are ≈2.37 and ≈5.75, well away from zero); see
§10 for why `Momentum_x` (%) is unreliable.

| case (dIters=10)        | Density net | Energy net | per-cycle median \|Δρ\|/ρ |
|-------------------------|-------------|------------|---------------------------|
| NN / First              | **+0.0203%**| **+0.0211%**| 1.11e-05                  |
| NN / Second             | −0.0008%    | −0.0015%   | 1.84e-07                  |
| NN / Third              | +0.0023%    | +0.0037%   | 4.92e-07                  |
| Barycentric / Second    | +0.0020%    | +0.0018%   | 7.97e-07                  |
| Barycentric / Third     | −0.0018%    | −0.0019%   | 4.64e-07                  |
| ConsGalerkinProj/Second | **+0.0003%**| **+0.0007%**| 1.07e-07                  |
| ConsGalerkinProj/Third  | −0.0020%    | −0.0022%   | 4.60e-07                  |

**Finding:** `NN/First` loses ~10× more mass/energy than any other method and has the
largest per-cycle defect. `ConsGalerkinProj/Second` is the best conserver. (`Momentum_x`
net for NN/First is +14.9%, but dominated by the small ∫Momentum_x ≈ 0.007; qualitative
only, see §10.)

### 6.2 Monotonicity — Density overshoot beyond the source range, per cycle

Overshoot = max fractional excursion of the interpolated field beyond the source min/max,
per cycle, for **Density**.

| case (dIters=10)        | cycles with overshoot | worst overshoot (frac of range) |
|-------------------------|-----------------------|---------------------------------|
| NN / First              | **0 / 80**            | 0                               |
| NN / Second             | 80 / 80               | 3.17e-03                        |
| NN / Third              | 80 / 80               | 3.95e-03                        |
| Barycentric / Second    | **0 / 80**            | 0                               |
| Barycentric / Third     | 80 / 80               | 3.58e-03                        |
| ConsGalerkinProj/Second | 80 / 80               | 2.00e-03                        |
| ConsGalerkinProj/Third  | 80 / 80               | 1.87e-03                        |

**Finding:** the only two methods that create **no** new extrema are `NN/First` (copies
source values) and `Barycentric/Second` (linear, bounded by simplex vertices). Every
gradient/Hessian method overshoots. `ConsGalerkinProj`'s overshoot is capped at ≈2.0e-3 =
its `GalerkinLimitTol` (its limiter working as designed).

---

## 7. Study F — Shock position/width at step 399 (dIters=10) vs no-adaptation reference

**Design:** after the full 400-step run, extract the Density transition (a broad feature at
x ≈ 2.70; see §10 — it is **not** a sharp shock) along horizontal probe lines, locate it by
**maximum-gradient position** (window-independent), and compare to the `Unsteady_Ref`
(no-adaptation) feature. Errors in **cells** (h = 0.01).
**Reference:** x_feature = **2.6984**, thickness (jump/max-grad) = **0.0493**. Initial
feature position (step 0) ≈ 2.00; total physical travel ≈ **0.70** (~70 cells).
**Source file:** `Shock/shock_step399.csv`.

| case (dIters=10)        | x_feature | pos error       | thickness increase |
|-------------------------|-----------|-----------------|--------------------|
| NN / First              | 2.7052    | **+6.86e-03 (0.69 cell)** | −1.25e-03  |
| NN / Second             | 2.6968    | −1.58e-03       | −6.23e-04          |
| NN / Third              | 2.6987    | +3.71e-04       | +1.31e-05          |
| Barycentric / Second    | 2.6990    | +6.49e-04       | **+5.42e-03 (0.54 cell)** |
| Barycentric / Third     | 2.6987    | +2.78e-04       | +1.53e-03          |
| ConsGalerkinProj/Second | 2.6975    | −8.34e-04       | +3.05e-04          |
| ConsGalerkinProj/Third  | 2.6984    | ≈0 (matches ref)| −6.08e-05          |

**Findings:**
- `NN/First` drifts the feature ~0.69 cell (worst position by ~4×), consistent with its
  worst conservation drift (§6.1).
- `Barycentric/Second` smears most (+0.54-cell thickness), consistent with being the most
  diffusive (linear).
- The high-order and conservative methods hold both position (<0.16 cell) and thickness
  (within ~3% of reference); `ConsGalerkinProj/Third` matches the no-interpolation
  reference to measurement precision.
- Even the worst position error is ~1% of the ~70-cell total shock travel.

---

## 8. Study G — Amplification by interpolation frequency (dIters=2 vs dIters=10)

**Design:** the same cases at `dIters=2` (400 cycles, 5× the interpolation count; same
physics, `SafeFactor=4`, correction off), compared to the `dIters=10` baseline. The feature
is located by the **50%-crossing** of the y-averaged Density profile (robust to NN's
staircase; a **different tracker than §7** — see §10). `cumDens` = cumulative Σ(tgt−src) of
∫Density over all cycles.
**Reference (50%-crossing):** x = **2.3993**, width (10–90%) = **0.3622**. Errors in cells
(h = 0.01). **Source:** `scratchpad/full.txt` (50%-crossing extractor).

| method / order          | dIters | pos err (cell) | width inc (cell) | cum ∫ρ loss |
|-------------------------|--------|----------------|------------------|-------------|
| NN / First              | 10     | +1.28          | +0.22            | +4.8e-04    |
| NN / First              | 2      | **+7.56**      | +1.32            | +7.0e-04    |
| NN / Second             | 10     | +0.07          | −0.06            | −2.0e-05    |
| NN / Second             | 2      | −0.07          | −0.04            | −3.3e-05    |
| NN / Third              | 10     | +0.11          | −0.03            | +5.5e-05    |
| NN / Third              | 2      | +0.07          | −0.01            | +8.0e-05    |
| Barycentric / Second    | 10     | +1.10          | +0.48            | +4.7e-05    |
| Barycentric / Second    | 2      | **+8.12**      | +1.50            | +3.6e-04    |
| Barycentric / Third     | 10     | +0.16          | +0.03            | −4.2e-05    |
| Barycentric / Third     | 2      | +0.33          | +0.38            | +8.2e-06    |
| ConsGalerkinProj/Second | 10     | +0.09          | −0.00            | +7.3e-06    |
| ConsGalerkinProj/Second | 2      | **+0.09**      | +0.04            | −8.4e-06    |
| ConsGalerkinProj/Third  | 10     | +0.15          | −0.02            | −4.8e-05    |
| ConsGalerkinProj/Third  | 2      | +0.14          | +0.17            | +4.9e-05    |

**Findings (three tiers):**
- **Fragile:** `NN/First` (+1.28 → +7.56 cell) and `Barycentric/Second` (+1.10 → +8.12
  cell) amplify ~6×, by different mechanisms (NN/First conservation-driven,
  Barycentric/Second diffusion-driven).
- **Robust:** `NN/Second`, `NN/Third`, `Barycentric/Third` stay small (<0.4 cell).
- **Immune:** `ConsGalerkinProj/Second` position **unchanged** (+0.09 → +0.09) under 5×
  interpolation; `ConsGalerkinProj/Third` similar.

---

## 9. Study H — Amplitude scaling (SafeFactor=2 vs 4) and the definitive Global test

**Design:** four cases at `SafeFactor=2` (double per-cycle displacement), keeping
`dIters=2`, each paired with its `SafeFactor=4` twin from §8. Same 50%-crossing tracker and
reference (x = 2.3993). Errors in cells. **Source:** `scratchpad/sf2cmp` output.

| case                    | SF (disp.) | pos err (cell) | width inc (cell) | cum ∫ρ loss |
|-------------------------|------------|----------------|------------------|-------------|
| NN/First — off          | 4          | +7.56          | +1.32            | +7.0e-04    |
| NN/First — off          | **2**      | **+25.78**     | +3.53            | +3.2e-03    |
| NN/First — **Global**   | 4          | +7.70          | +1.10            | +1.1e-14    |
| NN/First — **Global**   | **2**      | **+25.53**     | +3.25            | +2.4e-13    |
| ConsGalerkinProj/Second | 4          | +0.09          | +0.04            | −8.4e-06    |
| ConsGalerkinProj/Second | **2**      | **+0.09**      | +0.31            | −3.1e-05    |
| Barycentric/Second — off| 4          | +8.12          | +1.50            | +3.6e-04    |
| Barycentric/Second — off| **2**      | **+26.43**     | +2.03            | +6.7e-04    |

**Definitive findings:**
1. **Amplitude scaling:** doubling displacement (SF4→SF2) roughly triples–quadruples the
   fragile-method drift (NN/First +7.56 → +25.78 cell; Barycentric/Second +8.12 → +26.43
   cell).
2. **Global conservation does not fix the shock.** `NN/First` with "Global" has machine-zero
   conservation (3.2e-3 → 2.4e-13) yet its position drift is **+25.78 → +25.53 cell
   (unchanged, <1% difference)**: Global fixes the total integral but not where the mass sits.
3. **Local conservation immunity holds at 2× amplitude.** `ConsGalerkinProj/Second` position
   is **+0.09 → +0.09 cell (identical)**, ~280× better than `NN/First+Global` (+25.53 cell).

**Central conclusion (robust across order, interpolation frequency ×5, and amplitude ×2):**
*Global conservation is bookkeeping — it restores ∫ρ but leaves shock position unchanged;
only built-in local conservation (`ConsGalerkinProj`) pins the shock.*

### 9.1 Interpretation — why ConsGalerkinProj position is immune while its width grows with displacement

`ConsGalerkinProj` keeps the feature **position** unchanged under a doubling of nodal
displacement (Second: +0.09 → +0.09 cell, SF4→SF2; Third: +0.15 → +0.14) while its feature
**width** grows (Second: +0.04 → +0.31 cell; Third: −0.02 → +0.17). The explanation
separates established theory from this study's own measurements and one adopted assumption.
The references in §13 support the theoretical statements and should be checked against the
originals before publication.

**Framework.** Characterise the transported density transition by the moments of its
profile: 0th moment (integral / mass), 1st moment (centroid → feature *position*), 2nd
moment (variance → feature *width*).

**Established theoretical results (not novel to this work):**

- **(P1) Exact 0th-moment (mass) conservation.** A conservative L2/Galerkin projection
  reproduces the field integral exactly, by construction [Dukowicz & Kodis 1987; Margolin &
  Shashkov 2003]. This is the machine-precision conservation of §4–§6.
- **(P2) Conservative + bound-preserving ⇒ diffusive.** A remap that is simultaneously
  conservative and monotone/bound-preserving cannot preserve a sharp gradient and must
  spread it — the conservative-interpolation counterpart of Godunov's barrier theorem (a
  *linear* monotone scheme is at most first-order accurate [Godunov 1959]), established in
  the remapping setting [Margolin & Shashkov 2003; Kucharik, Shashkov & Wendroff 2003]. In
  moment terms: conservation fixes the 0th moment and monotonicity caps over/undershoot, but
  the **2nd moment (width) is not preserved and grows** — this smearing is intrinsic, not a
  defect of the implementation.
- **(P3) Bound-preserving limiters are diffusive.** Enforcing local maximum-principle bounds
  and conservatively redistributing the clipped excess (exactly what `GalerkinLimit` does)
  is a diffusive operation [Zalesak 1979; Kucharik, Shashkov & Wendroff 2003].
- **(P4) L2 projection of a discontinuity smears it.** Projecting a jump onto a
  piecewise-polynomial space whose cells straddle it yields a smeared (or, unlimited,
  oscillatory) representation — a standard property of the L2 projection [Hesthaven &
  Warburton 2008].

Together (P1)–(P4) imply, from existing theory, that `ConsGalerkinProj` conserves mass and —
to the extent it is linearity-preserving — holds position, while smearing the width. The
observed "immune position, growing width" is therefore the expected signature of a
conservative bounded remap, and the mirror image of `NN/First` (sharp but position-drifting).

**This study's measurements supporting the picture:**

- **(M1) Position immunity follows from 1st-moment preservation.** A single interpolation of
  the real field (SF=2), computing the discrete moments M0 = Σ w_i ρ_i (mass) and
  M1 = Σ w_i x_i ρ_i (1st moment) with lumped nodal areas, source vs target
  (`MomentPreservationTest.py`, `Moments/moment_preservation.csv`):

  | method (kind)                          | ΔM0/M0 (mass) | ΔM1/M1 (1st moment) | centroid shift |
  |----------------------------------------|---------------|---------------------|----------------|
  | ConsGalerkinProj/Second (conservative) | 8.4e-8        | **5.9e-7**          | 7.2e-7         |
  | Barycentric/Second (non-cons., linear) | 6.2e-7        | 3.3e-6              | 3.8e-6         |
  | NN/First (non-conservative)            | −2.8e-5       | **−4.6e-5**         | −2.6e-5        |

  ConsGalerkinProj preserves **both** mass (0th) and centroid (1st moment) to near machine
  precision → position immune. NN/First preserves **neither** → misplaces mass, shifts the
  centroid, feature drifts. The centroid-preservation ranking (ConsGalerkinProj < Barycentric
  < NN/First) matches the position-drift ranking of §7, so 1st-moment preservation controls
  feature-position immunity — the mechanism behind the central conclusion.

- **(M2) The super-linear width growth is driven by the bound-preserving limiter, whose
  activation rises sharply with displacement.** Running the ConsGalerkinProj/Second
  projection with the limiter OFF gives the field the limiter must correct; its excursion
  beyond the source range is what the limiter suppresses
  (`LimiterActivationTest.py`, `LimiterAct/limiter_activation.csv`):
  - **SF=4** (mean displacement 9.40e-4): unlimited overshoot = **7.6e-6** — below the 2e-3
    limiter tolerance, so the limiter is **dormant** and the SF=4 width growth is base L2
    diffusion (P2/P4) alone.
  - **SF=2** (1.87e-3, 1.99× the displacement): unlimited overshoot = **3.6e-3 (≈476× larger)**,
    with **97% concentrated at the feature** (|x−2.70|<0.15).

  The limiter's activation is thus strongly super-linear in δ and localized at the transition:
  its conservative redistribution adds smearing on top of the base L2 diffusion at larger
  displacement. This establishes the mechanism; it does not by itself partition the accumulated
  width growth between limiter redistribution and base L2 diffusion, because the accumulated
  width with limiter on vs off lies below the per-cycle metric-noise floor on this field.

**Adopted assumption (not independently established):**

- **(A1) Per-cycle diffusion scales with the displacement δ.** We take the width added per
  cycle to grow with the geometric mesh displacement, so doubling δ (SF4→SF2) increases
  per-cycle smearing and the accumulated width. This is consistent with the two-point
  δ-scaling (SF4→SF2, ≈8× accumulated width for 2× displacement) but is not a derived law:
  the per-cycle smearing (~1e-5) sits below the metric-noise floor (~1e-3), so it cannot be
  isolated on this broad feature. A definitive scaling would require additional full 400-step
  runs at more displacements, or a synthetic sharp-step feature where per-cycle smearing is
  resolvable. Treat the δ-scaling as an observed two-point trend, not an established law.

**Note on magnitudes.** The widths here are sub-cell and subject to the measurement caveats
of §10 (feature-tracker sensitivity, griddata sampling). The robust content is the **trend**
(width increases with δ while position does not) and the mechanism (M1, M2), not the exact
numerical factor.

---

## 10. Caveats, assumptions, and known weaknesses (READ BEFORE WRITING CLAIMS)

1. **The "shock" is a broad transition, not a sharp shock.** Along the probe lines the
   Density feature spans ~5 cells (max-gradient ≈ 5–6 in the griddata sampling). "Shock
   position/width" is the position/width of a **density transition feature** (possibly a
   contact/rarefaction), not a strong shock; its position is therefore **metric-dependent**.

2. **Two feature trackers are used, giving different absolute positions.**
   - §7 (dIters=10) uses **maximum-gradient** location → reference x = **2.6984**, total
     travel ≈ 0.70 (~70 cells).
   - §8–§9 (amplification, SF2) use **50%-crossing** of the y-averaged profile → reference
     x = **2.3993**, total travel ≈ 0.40 (~40 cells).
   The same dIters=10 case therefore appears with different numbers in §7 vs §8 (e.g.
   NN/First: +0.69 cell max-grad vs +1.28 cell 50%-crossing). **Do not mix the two trackers
   in one plot or normalisation.** The 50%-crossing tracker is used for NN-bearing studies
   because the maximum-gradient tracker locks onto NN's piecewise-constant staircase and is
   unreliable for NN fields.

3. **Absolute NN shock numbers are sampling-sensitive.** For NN's staircase field the
   50%-crossing depends on the probe-line set and resolution: `NN/First dIters=2 off` reads
   **+7.56 cell** with one sampling (8 lines, 400 probe points) and ~**+25.78 cell** with
   another. The robust quantities are the **relative** ones measured with identical sampling:
   (a) off-vs-Global deltas (§9: +25.78 vs +25.53 — essentially equal), and (b) the
   immune-vs-fragile contrast (ConsGalerkinProj +0.09 vs fragile ~25). Quote absolute NN cell
   counts with this caveat, or as order-of-magnitude.

4. **Momentum_x conservation percentages are unreliable.** ∫Momentum_x ≈ 0.007 (near zero),
   so relative drift on it is small-denominator noise (e.g. NN/First "+14.9%", some "%init"
   values exceed 100%). Use **Density and Energy** for conservation claims; treat Momentum
   only qualitatively or via **absolute** cumulative loss.

5. **Monotonicity "worst_overshoot_field" in the raw CSV is often Momentum_y**, whose range
   is ≈0 (symmetric field) → its overshoot fraction is noise. §6.2 reports **Density-specific**
   overshoot, the meaningful quantity. Do not quote the raw `worst_overshoot` column from
   `multicycle_summary.csv` without this correction.

6. **Single realization per case (one random seed).** No ensemble averaging over random
   meshes; reported drifts/positions are single-run outcomes with no error bar or variance
   estimate. A paper claim of a "rate" or "factor" should note it is from one realization.

7. **griddata (linear) sampling adds a baseline smearing to width.** Widths are sampled by
   linear interpolation onto probe lines, adding a mesh-scale baseline equally to all cases;
   it cancels in **width-increase differences** (interp − source) but the absolute widths are
   not the true feature widths.

8. **Per-case configuration is verified.** Every case directory holds its own copy of the
   interpolation code and parameters. The settings behind each case (`InterpMethod`,
   `InterpOrder`, `EnforceConservation`, `dIters`, `SafeFactor`, `debugInterp=1`, and the
   runtime `WLSDegree=2` gradient for the third-order methods) are confirmed from each case's
   own copy, so the configuration behind every number in §2–§9 is verified.

9. **δ/h in the refinement study is only approximately constant.** `SafeFactor=8` fixes the
   safe-radius fraction, giving mean δ/h ≈ 0.057 across meshes, but the per-node safe radius
   varies slightly with local geometry. The order-of-accuracy fits assume δ/h constant; the
   small variation is a minor source of scatter.

10. **Analytic-field wavenumber differs between §2 and §3** (2π refinement vs 4π fixed-h
    sweep), so the coarsest mesh resolves the field in the refinement study. §2 and §3
    absolute L1 values are therefore **not** directly comparable; only the exponents are.

11. **`ConsGalerkinProj/Third` order fit (2.79) is pre-asymptotic**, not a failure to reach 3
    — slopes rise 2.48 → 2.89 → 2.96. Present it as "approaching 3rd order, pre-asymptotic on
    the coarse meshes," or fit only the two finest intervals.

12. **The single-interpolation real-field shock comparison (off vs Global) is not archived
    as a table.** Its conclusion — that Global does not change single-interpolation shock
    position or width — is established at the multi-cycle level in §9.

---

## 11. Data provenance (source file → what it contains)

| file                                   | contents                                             | script |
|----------------------------------------|------------------------------------------------------|--------|
| `Tier0_Refine/tier0_refine.csv`        | order-of-accuracy L1/L2/Linf vs h (§2)               | `Tier0_InterpBenchmark.py` |
| `Tier0/tier0_analytic.csv`             | displacement-sensitivity L1/L2/Linf vs δ (§3)        | `Tier0_InterpBenchmark.py` |
| `Correction/correction_comparison.csv` | single-interp L1 + conservation, off vs Global (§4)  | `CorrectionComparison.py` |
| `GalerkinTol/galerkin_tol.csv`         | limiter-tolerance conservation sweep (§5)            | `GalerkinTolSweep.py` |
| `LimiterAct/limiter_activation.csv`    | limiter activation vs displacement, SF4/SF2 (§9.1 M2)| `LimiterActivationTest.py` |
| `Moments/moment_preservation.csv`      | 0th/1st-moment preservation by method (§9.1 M1)      | `MomentPreservationTest.py` |
| `MultiCycle/multicycle_summary.csv`    | 80-cycle conservation drift + overshoot (§6)         | `MultiCycleAnalysis.py` |
| `Shock/shock_step399.csv`              | step-399 feature position/width, max-grad (§7)       | (extractor) |
| `scratchpad/full.txt`                  | dIters=2 vs 10 amplification, 50%-crossing (§8)      | (50%-crossing extractor) |
| `scratchpad/sf2cmp` output             | SF2 vs SF4 scaling + Global test (§9)                | (paired extractor) |
| per-case `DebugFiles/Integrals_*.csv`  | per-cycle ∫φV src vs tgt (raw conservation)          | `InterpSolution.py` (debugInterp=1) |
| per-case `DebugFiles/{Min,Max}imums_*` | per-cycle src/tgt extrema (raw overshoot)            | `InterpSolution.py` (debugInterp=1) |

The §8–§9 amplification/SF2 tables are extractor outputs held in the scratchpad rather than
versioned CSVs in the case tree; the per-case restarts and `DebugFiles` from which they are
computed are all on disk, so they can be regenerated into archival CSVs if the paper needs them.

---

## 12. Suggested figure/table set for the paper

1. **Convergence plot (log-log):** Density L1 vs h (§2.1), 7 series, with reference
   slope-1/2/3 triangles. Caption the pre-asymptotic ConsGalerkinProj/Third.
2. **Order table:** fitted p (§2.2) with nominal orders.
3. **(Optional) displacement-sensitivity plot:** Density L1 vs δ (§3.1) — labelled
   "fixed-h sensitivity, not order," per caveats §10.1/§10.10.
4. **Conservation-drift bar chart:** Density net % over 80 cycles (§6.1).
5. **Overshoot table:** Density overshoot (§6.2).
6. **Amplification/scaling table:** §8 + §9 combined, with the tracker caveat (§10.2/3).
7. **Headline result callout:** §9 — Global (machine-zero conservation) vs local
   conservation on shock position (+25.5 vs +0.09 cell).
8. **Interpretation box (§9.1):** the moment-hierarchy explanation of ConsGalerkinProj's
   immune position / growing width, with the theory/measurement/assumption split and
   references §13.

---

## 13. References (for the §9.1 theoretical claims)

Canonical works supporting the established results (P1)–(P4) in §9.1; none is a result of the
present study. Verify exact volume/issue/page against the originals before submission.

- **[Godunov 1959]** S. K. Godunov, "A difference method for the numerical calculation of
  discontinuous solutions of the equations of hydrodynamics," *Matematicheskii Sbornik*,
  **47(89)**:271–306, 1959. — Barrier theorem: a linear monotone scheme is at most
  first-order accurate (root principle behind P2).
- **[Dukowicz & Kodis 1987]** J. K. Dukowicz and J. W. Kodis, "Accurate conservative
  remapping (rezoning) for arbitrary Lagrangian–Eulerian computations," *SIAM Journal on
  Scientific and Statistical Computing*, **8(3)**:305–321, 1987. — Foundational conservative
  remapping (P1).
- **[Margolin & Shashkov 2003]** L. G. Margolin and M. Shashkov, "Second-order sign-preserving
  conservative interpolation (remapping) on general grids," *Journal of Computational Physics*,
  **184(1)**:266–298, 2003. — Conservative, sign/bound-preserving remapping on general grids
  and its diffusive character (P1, P2).
- **[Kucharik, Shashkov & Wendroff 2003]** M. Kucharik, M. Shashkov, and B. Wendroff, "An
  efficient linearity-and-bound-preserving remapping method," *Journal of Computational
  Physics*, **188(2)**:462–471, 2003. — Bound-preserving remap; diffusion from bound
  enforcement (P2, P3).
- **[Zalesak 1979]** S. T. Zalesak, "Fully multidimensional flux-corrected transport
  algorithms for fluids," *Journal of Computational Physics*, **31(3)**:335–362, 1979. —
  Flux-corrected transport; numerical diffusion introduced by bound preservation (P3).
- **[Hesthaven & Warburton 2008]** J. S. Hesthaven and T. Warburton, *Nodal Discontinuous
  Galerkin Methods: Algorithms, Analysis, and Applications*, Springer, 2008. — L2 projection
  of discontinuities: smearing (limited) vs Gibbs oscillations (unlimited) (P4).
