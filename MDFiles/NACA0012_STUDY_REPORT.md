# NACA0012: Mesh-to-Mesh Interpolation and BDF2 Restart

## Summary

A steady, shock-free flow is used as a controlled setting for two questions about
transferring a solution between meshes. Because the flow is inherently steady, any
departure of the force coefficients from an independently computed reference is
attributable to the transfer.

**Interpolation error does not average out for drag.** Under repeated interpolation onto
a randomly re-moved mesh, the drag trend reproduces across independent realisations of
the random field to within 7%, while the lift trend does not reproduce at all. Lift and
moment errors are zero-mean and cancel over many cycles; interpolation dissipation is
sign-definite, and drag is the coefficient that integrates it. Local conservation reduces
the drag drift sixfold without removing it.

**The perturbation reaches the airfoil as an acoustic wave.** A single interpolation
confined to a narrow annulus produces a force response whose delay is linear in the
launch radius, with zero intercept. Probing the pressure field directly gives a
propagation speed of 1.18 ± 0.08 chords per time-unit against a sound speed of 1.18322.
Doubling the perturbation amplitude at fixed radius moves the arrival by at most two
timesteps, so the delay is propagation and not a detection threshold.

**A BDF2 restart cannot be characterised by the error in its latest state.** A defect
placed only in the older of the two history states has zero latest-state error, yet moves
the lift coefficient by ten to thirteen times the native noise floor, with the sign the
−1/2 history coefficient requires. The measured initial-residual difference matches the
predicted history forcing to export precision for every transfer method.

**Global conservation enforcement does not change the dynamic response.** It drives the
density integral defect to a common value while leaving the state defect, the first-step
response and the lift response unchanged in the fourth significant figure.

**A silent defect in the quadratic projection is corrected here.** Its supermesh
double-counted the overlap of seven cells, producing a plausible restart with no warning.
A constant-field patch test exposes it immediately, and with the correction the
third-order Galerkin transfer becomes the most accurate of the twelve.

---

## 1. Case

Inviscid Euler flow over a NACA0012 at **M = 0.1**, **AoA = 6°**. Deeply subsonic, so the
field is smooth and shock-free; this is deliberately not a shock-capturing test.
Non-dimensionalisation sets p = ρ = T = 1, giving a sound speed of

    a = sqrt(1.4) = 1.18322 chords per time-unit

and a free-stream speed of 0.118322. Time integration is second-order dual-time stepping
at **dt = 1e-3**, so acoustic travel is 1.18e-3 chords per timestep and a perturbation
released at radius d needs

    lag = d / (a·dt) = 845·d timesteps

to reach the surface. Convection is ten times slower, so the two timescales are cleanly
separated. Forces use `MARKER_MONITORING = (airfoil)`, REF_AREA = REF_LENGTH = 1, and a
moment reference at x = 0.25; the flow being inviscid, they are pressure-only.

## 2. Meshes

Two mesh pairs are used, each suited to its question.

**Part I** perturbs a single mesh. `Coarse_Tria_Euler.su2` — 16 255 nodes, 32 077
triangles, farfield at 100 chords — has nodes displaced inside a chosen annulus measured
by wall distance, with the boundary held fixed. The two meshes are then identical outside
the annulus, which localises the transfer disturbance to the launch band exactly.

**Part II** uses two independently generated meshes.

| | A: `Coarse_Tria_Euler.su2` | B: `Coarse_Tria_Euler_Delaunay.su2` |
|---|---|---|
| nodes / triangles | 16 255 / 32 077 | 22 508 / 44 583 |
| farfield | 100 chords | 100 chords |
| first off-wall spacing | 5.9e-4, first-layer median 3.7e-3 | 4.3e-4, first-layer median 2.7e-3 |
| nodes within 0.01 chords of the wall | 5.2% | 5.6% |
| edge-ratio median / 99th percentile | — | 1.19 / 1.57 |

Mesh B has no inverted, duplicate or orphan elements, and its total area matches the
analytic area of the 63-sided farfield polygon less the airfoil section to 1.7e-8
relative. The volume discretisations are independent: only the boundary nodes coincide,
and the median distance from a node of B to the nearest node of A is 1.6e-2.

The pair shares its boundary discretisation — the same 370 airfoil nodes and 63-node
farfield polygon — so the discretised body is identical and ΔC_p is a node-to-node
comparison on common surface points. The defect nonetheless has full surface support: at
a shared boundary node the transfer is exact, so the defect there is the native
source/target difference, which reaches 3.5e-2.

## 3. Conventions

**Mass matrix.** M is diagonal and holds the median-dual control volume of each node. On
a simplex the median dual splits every triangle into three equal-area regions, so a
node's volume is one third of the area of each incident triangle. Summed volumes close on
the total mesh area to 3.6e-12.

**Residual.** The solver stores the negative of the residual used here,

    LinSysRes = −Φ(U^{n+1}; U^n, U^{n-1}),
    Φ = (1/Δt)·M·(3/2·U^{n+1} − 2·U^n + 1/2·U^{n-1}) + F(U^{n+1}),

and the exported field is the complete BDF2 residual, time term included, evaluated at
the trial state before the first update. Setting U^n = U^{n-1} = U makes the time term
vanish identically, which is how the spatial residual F(U) is recovered.

**Transfer methods.** `NN` is a nearest-neighbour point-cloud method with a Taylor
correction from weighted-least-squares gradients; `Barycentric` locates the containing
simplex; `ConsGalerkinProj` is an L2 projection over the exact supermesh, locally
conservative by construction. Orders First, Second and Third select the reconstruction
degree. `EnforceConservation = "Global"` rescales the result to restore chosen total
integrals.

---

# Part I — Interpolation onto a randomly re-moved mesh

## 4. Force drift under repeated interpolation

Interpolation is applied every 10 timesteps over a 400-step run — 39 events — with the
perturbation annulus spanning 0.01 to 20 chords, so essentially the whole domain moves.
The reference is the same case run with no adaptation and no interpolation, which holds
the coefficients flat. Two methods are compared at second order, each with two
independent seeds of the random displacement field.

Trends are least-squares slopes of the deviation from the reference, quoted per 100
iterations with the standard error of the slope, so "no drift" means the slope is not
resolved from zero rather than merely small.

| | ConsGalerkinProj | | NN | |
|---|---|---|---|---|
| | seed 1 | seed 2 | seed 1 | seed 2 |
| CL trend /100it | −2.4e−6 ± 6.1e−6 | +2.1e−6 ± 5.7e−6 | +1.2e−5 ± 7.0e−6 | +4.9e−5 ± 8.4e−6 |
| CD trend /100it | **+6.8e−6 ± 2.3e−6** | **+6.4e−6 ± 2.2e−6** | **+3.9e−5 ± 3.3e−6** | **+3.8e−5 ± 3.5e−6** |
| CMz trend /100it | +1.3e−7 ± 1.2e−6 | −1.0e−6 ± 1.1e−6 | **−1.46e−5 ± 1.5e−6** | **−1.49e−5 ± 1.7e−6** |
| CD scatter | 5.3e−5 | 5.1e−5 | 8.8e−5 | 9.2e−5 |

Reproducibility across seeds is what separates a systematic drift from a realisation. The
drag trend reproduces to 7% for the conservative method and 3% for the non-conservative
one, and the moment trend of the non-conservative method to 2%. The lift trend does not
reproduce: it changes sign between seeds for `ConsGalerkinProj` and by a factor of four
for `NN`. A linear fit through a bounded random walk returns a non-zero slope on any
single realisation, which is what the lift signal is.

In relative terms the drag climbs **+0.065 %/100it** under `ConsGalerkinProj` and
**+0.389 %/100it** under `NN`, reaching roughly +0.6% and +1.8% of CD by step 400.

The asymmetry follows from the sign structure of the errors. Lift and moment respond to
interpolation error that is as likely positive as negative, and random re-meshing averages
it away over 39 cycles. Interpolation dissipation is not sign-symmetric: each projection
smooths the field, generating entropy that never comes back, and drag is the surface
integral most sensitive to it. This also explains why global conservation enforcement is
not the remedy — `ConsGalerkinProj` already conserves the integral exactly, and its drag
still climbs.

No claim is made about the time-dependence of the trend within the run; only the mean
slope over the full 400 steps reproduces across seeds.

## 5. Delay scaling with perturbation distance

A single interpolation is fired at one step into a narrow annulus at radius d, with the
run continuing afterwards with no further adaptation. The loop structure ties the response
window to the event step, so a window of W steps requires the event at step W. Because
every onset is referred to the case's own peak, the event step is chosen greater than
twice the predicted lag, leaving room for the peak at roughly 1.2·845·d and for the pulse
to decay behind it: step 200 for d = 0.05 and 0.10, step 400 for d = 0.20.

Up to the event the perturbed run is **bit-identical** to the reference — maximum
deviation 0.0e+00 in CL and CD over the preceding 198 steps — so the measurement floor is
machine precision and the signal sits five orders of magnitude above it.

The response is a damped oscillating pulse, so onset is detected as the first crossing of
a fixed fraction of that case's own peak, which normalises out the amplitude differences
between bands. Several thresholds are reported because the choice is arbitrary.

| d | onset 1% | 5% | 25% | 50% | peak lag | predicted 845·d |
|---|---|---|---|---|---|---|
| 0.05 | 16 | 19 | 26 | 31 | 74 | 42 |
| 0.10 | 37 | 44 | 75 | 82 | 97 | 85 |
| 0.20 | 105 | 115 | 144 | 171 | 203 | 169 |

Fitting lag = C·d + b:

| threshold | C | intercept | R² | C / 845 |
|---|---|---|---|---|
| 1% | 606 | −18.0 | 0.989 | 0.72 |
| 5% | 650 | −16.5 | 0.994 | 0.77 |
| 25% | 773 | −8.5 | 0.992 | 0.91 |
| 50% | 927 | −13.5 | 0.999 | 1.10 |

The functional form is the result, not the coefficient. A wave gives lag proportional to d
with zero intercept; diffusion would give d²; an instantaneous response would give a lag
independent of d. The measurement is linear with a near-zero intercept at every threshold,
and the acoustic prediction is bracketed between the 25% and 50% detectors. A
fraction-of-peak crossing on a smoothly rising ramp fires before the front, so low
thresholds return a lower bound on the lag.

**Amplitude control.** Doubling the displacement at fixed launch radius raises the peak
response by a factor 1.63 and shifts the arrival by **+1, −1, −1 and −2 timesteps** at the
1%, 5%, 25% and 50% thresholds. A delay that were really a detection threshold would
shrink markedly under a 63% stronger pulse. The delay is propagation.

## 6. Propagation speed

The force is a surface integral and superposes contributions arriving from the whole
annulus over a spread of travel times, which smears the front. Probing the pressure field
removes that limitation.

Everything inside r < d is undisturbed at the moment of the event and its nodes never
move, so a fixed node index is a fixed point in space. A pulse launched at radius d
reaches a shell at radius r after 845·(d − r) timesteps, and sweeping r inside one case
traces the propagation directly. A shell is equidistant from the launch annulus all the
way round, so its arrival is sharp. The fit is differential — regressing lag on (d − r)
uses only differences between shells, so any constant bias in the onset detector falls
into the intercept and leaves the slope untouched. The measured quantity is the shell-RMS
of the pressure difference, which averages azimuthally.

Shells within 0.08 chords of the annulus respond at lag 0 to 2 regardless of radius: they
lie inside the interpolation's own error footprint and are perturbed directly rather than
reached by a wave. They are excluded from the fit.

| threshold | S | R² | implied speed |
|---|---|---|---|
| 2% | 907 | 0.960 | 1.103 |
| 5% | 862 | 0.975 | 1.160 |
| 10% | 917 | 0.993 | 1.091 |
| 25% | 826 | 0.999 | 1.210 |
| 50% | 794 | 0.998 | 1.259 |
| 75% | 782 | 0.998 | 1.279 |

Mean **S = 848 ± 57 against the predicted 845**, a propagation speed of **1.18 ± 0.08
against the sound speed 1.18322**. The six estimates span 782 to 917, so the speed is
determined to roughly ±7%, and no threshold choice pulls it away from the acoustic value.

**Choice of reference.** The measurement above is referenced to the unperturbed
trajectory on the source mesh. That is the correct baseline here, and not a confound,
because the two meshes are identical outside the annulus: an L2 projection between
identical meshes is the identity, so the transfer is exact there and the disturbance is
genuinely confined to the launch band. Referencing instead to an independently converged
trajectory on the deformed mesh introduces the native source/target discretisation
difference — 2.3e-5 in CL — as a static offset present at the airfoil from the instant of
transfer, which is not interpolation error. That offset saturates the lower detection
thresholds; the propagation is recovered at the unsaturated ones, giving S = 792 and 781
at the 50% and 75% thresholds with R² = 0.999, a speed of 1.26 to 1.28. Both references
therefore return the acoustic speed within the precision of the measurement.

**Causality under an implicit solver.** Dual-time stepping with a global linear solve at
every inner iteration has a numerical domain of dependence covering the whole mesh, so the
deviation is non-zero everywhere immediately and strict acoustic causality is not enforced
numerically. What survives is that this instantaneous component is tiny and collapses with
distance: the deviation at the event step itself is 6.4e-8 for d = 0.05 against 1.2e-9 for
d = 0.10 and d = 0.20, with the latter two at the output precision floor. Causality is
respected in amplitude rather than in support, which is why the wavefront cannot be
resolved by a first-non-zero criterion and why the differential slope is the right
statistic.

---

# Part II — BDF2 restart validation

## 7. Native solutions

Both steady solutions are converged to the residual floor of the shared configuration:
mesh A reaches rms[Rho] = −12.85 in 4402 iterations, mesh B reaches −13.00 in 3483. Each
is then advanced from its own steady state with both BDF2 histories equal, giving the
native trajectories against which every response is measured.

| | A | B | B − A |
|---|---|---|---|
| CL | 0.69908195 | 0.70902504 | +9.94e-3 |
| CD | 0.00997182 | 0.00305620 | −6.92e-3 |
| CMz | 0.00965569 | 0.00722540 | −2.43e-3 |

Coefficient variation over the 51-step native window on mesh B — the noise floor against
which every response is judged — is 1.9e-4 in CL, 1.4e-4 in CD and 5.7e-5 in CMz. The
window is a smooth relaxation transient rather than random scatter, and it cancels because
responses are differenced against the native trajectory step by step.

Inviscid drag is entirely spurious, so CD is physically zero and the 69% relative gap
between the meshes carries no meaning. Drag responses are reported as absolute changes.

## 8. History configurations and transfer defects

Second-order dual time consumes restarts at RESTART_ITER−1 and RESTART_ITER−2, so a
history pair is built by writing two restart files. With ε = Î_AB U_A* − U_B* and
η = 2ε^n − ε^{n-1}/2:

| case | U^n | U^{n-1} | ε^n | ε^{n-1} | η |
|---|---|---|---|---|---|
| `native` | U_B* | U_B* | 0 | 0 | 0 |
| `full` | Û_B* | Û_B* | ε | ε | 3ε/2 |
| `latest_only` | Û_B* | U_B* | ε | 0 | 2ε |
| `older_only` | U_B* | Û_B* | 0 | ε | −ε/2 |

| method | ‖ε‖_M | relative | integral defect, density |
|---|---|---|---|
| NN1 | 1.474e-3 | 3.08e-6 | −7.284e-4 |
| NN2 | 7.243e-4 | 1.51e-6 | −7.937e-4 |
| NN3 | 7.157e-4 | 1.50e-6 | −7.983e-4 |
| BAR2 | 7.040e-4 | 1.47e-6 | −8.136e-4 |
| BAR3 | 7.117e-4 | 1.49e-6 | −7.911e-4 |
| GP2 | 7.137e-4 | 1.49e-6 | −8.121e-4 |
| GP3 | 7.126e-4 | 1.49e-6 | −7.853e-4 |

with ‖U_B*‖_M = 478.46. Every method except NN1 lands within 3% of the same defect norm,
which is the native source/target discretisation difference; the transfer operator
contributes little beyond it. First-order nearest-neighbour doubles that norm. The
third-order Galerkin projection is the most accurate of the set, by a small margin
concentrated in the far field, and reaching that behaviour required the correction
described in §14.

## 9. Initial-residual identity

A physical step begins from U^n, so `full` and `latest_only` start at Û_B* while `native`
and `older_only` start at U_B*. Residual differences are therefore formed only within a
pair sharing U^n. Two such pairs exist and both are used: (`older`, `native`) at trial
U_B* with Δη = −ε/2, and (`latest`, `full`) at trial Û_B* with Δη = +ε/2 — the same
magnitude with opposite sign. Converged first-step responses are unaffected, since a
tightly converged step reaches the same fixed point regardless of its starting iterate.

| | relative identity error | cosine |
|---|---|---|
| the ten point-location methods | 2.4e-8 to 4.2e-8 | 1.0000012 to 1.0000017 |
| GP2 | 4.4e-6 | 1.0000016 |
| GP3 | 1.3e-6 | 1.0000008 |

Snapshot fields are stored in single precision, whose epsilon is 1.19e-7; differencing two
vectors of similar magnitude amplifies that, so the largest errors observed are at export
precision. The identity is satisfied for every method on both pairs.

Because the exported residual already contains the time term, this identity is an
implementation check on the solver's BDF2 assembly together with the restart writer and
the volume computation. It validates the history algebra, not the state response.

## 10. First-step response and the two-history argument

| method | case | ‖δ^{n+1}‖_M | ΔCL | ΔCD | ΔCMz |
|---|---|---|---|---|---|
| NN1 | full | 1.461e-3 | +5.765e-3 | +1.125e-2 | +1.307e-4 |
| | latest | 1.948e-3 | +7.683e-3 | +1.502e-2 | +1.737e-4 |
| | older | 4.872e-4 | −1.832e-3 | −3.629e-3 | −5.049e-5 |
| NN2 | full | 7.230e-4 | −6.907e-3 | +5.107e-3 | +1.877e-3 |
| | latest | 9.638e-4 | −9.198e-3 | +6.810e-3 | +2.495e-3 |
| | older | 2.415e-4 | +2.321e-3 | −1.710e-3 | −6.261e-4 |
| BAR2 | full | 7.032e-4 | −7.308e-3 | +6.297e-3 | +1.671e-3 |
| | latest | 9.374e-4 | −9.740e-3 | +8.392e-3 | +2.213e-3 |
| | older | 2.349e-4 | +2.455e-3 | −2.107e-3 | −5.580e-4 |
| GP2 | full | 7.128e-4 | −7.598e-3 | +6.134e-3 | +1.785e-3 |
| | latest | 9.502e-4 | −1.012e-2 | +8.173e-3 | +2.368e-3 |
| | older | 2.381e-4 | +2.547e-3 | −2.055e-3 | −5.986e-4 |
| GP3 | full | 7.116e-4 | −7.130e-3 | +5.471e-3 | +1.670e-3 |
| | latest | 9.486e-4 | −9.501e-3 | +7.292e-3 | +2.214e-3 |
| | older | 2.377e-4 | +2.395e-3 | −1.832e-3 | −5.566e-4 |

The `older_only` case is the decisive one. Its latest-state defect is exactly zero, so any
restart metric built on ε^n alone predicts no response. The measured lift response is
1.8e-3 to 2.5e-3 against a native noise floor of 1.9e-4 — ten to thirteen times the floor —
and its sign is opposite to the `full` case throughout, as the −1/2 coefficient on the
older history requires.

Surface diagnostics place the largest discrepancies at the leading and trailing edges: the
maximum |ΔC_p| for the `full` configuration is 0.45 at x = 0.001 for NN1, 0.40 at
x = 0.999 for BAR2, 0.29 at x = 0.003 for GP2 and 0.30 at x = 0.999 for GP3.

## 11. Superposition and the linear regime

| method | state superposition defect | latest/full | older/full | force superposition defect, CL | ΔCL latest/full | ΔCL older/full |
|---|---|---|---|---|---|---|
| NN1 | 1.81e-2 | 1.3333 | 0.3335 | 1.50e-2 | 1.3328 | −0.3178 |
| NN2 | 3.52e-2 | 1.3331 | 0.3340 | 4.21e-3 | 1.3318 | −0.3360 |
| NN3 | 3.40e-2 | 1.3331 | 0.3340 | 5.29e-3 | 1.3309 | −0.3362 |
| BAR2 | 3.40e-2 | 1.3332 | 0.3341 | 3.06e-3 | 1.3328 | −0.3359 |
| BAR3 | 3.35e-2 | 1.3332 | 0.3340 | 2.84e-3 | 1.3329 | −0.3358 |
| GP2 | 3.40e-2 | 1.3330 | 0.3340 | 2.99e-3 | 1.3322 | −0.3352 |
| GP3 | 3.45e-2 | 1.3331 | 0.3341 | 3.27e-3 | 1.3327 | −0.3359 |

predicted 4/3 = 1.3333 and −1/3 = −0.3333.

The state ratios recover the predicted values to within 0.3% for every method, and the
force ratios to within 1% except for NN1, which is 4.7% low on the `older`/`full` lift
ratio. Superposition defects are 1.8% to 3.5% in the state and 0.3% to 1.5% in lift. The
transfer defect already lies inside the linear regime, so no amplitude sweep is required
to support the linearised interpretation.

## 12. Global conservation correction

Applied to the nearest-neighbour and barycentric transfers, with identical history
configurations and solver settings.

| pair | ‖ε‖_M | integral defect, density | ‖δ_full‖_M | ΔCL, full |
|---|---|---|---|---|
| NN1 | 1.47423e-3 | −7.28431e-4 | 1.46104e-3 | +5.76465e-3 |
| NN1 + Global | 1.47407e-3 | −8.12114e-4 | 1.46091e-3 | +5.76976e-3 |
| BAR2 | 7.03957e-4 | −8.13552e-4 | 7.03146e-4 | −7.30790e-3 |
| BAR2 + Global | 7.03956e-4 | −8.12114e-4 | 7.03146e-4 | −7.30770e-3 |

The correction drives the density integral defect to −8.12114e-4 for every method, the
value set by the native source/target integral difference, which is what it is designed to
do. Everything dynamic is untouched: the defect norm, the first-step response norm and the
lift response all agree with the uncorrected transfer in the fourth significant figure,
and the superposition defects and response ratios are unchanged. Removing the scalar
global defect does not remove the spatial defect that drives the restart.

## 13. Relaxation over the response window

Lift defect over the 50 steps following the restart.

| method | case | ΔCL first step | max abs, at step | ΔCL last step |
|---|---|---|---|---|
| NN1 | full | +5.765e-3 | 2.697e-2 at 5 | −1.390e-4 |
| NN2 | full | −6.907e-3 | 6.907e-3 at 2 | −2.702e-3 |
| BAR2 | full | −7.308e-3 | 7.308e-3 at 2 | −1.627e-3 |
| GP2 | full | −7.598e-3 | 7.598e-3 at 2 | −1.946e-3 |
| GP3 | full | −7.130e-3 | 7.130e-3 at 2 | −2.399e-3 |

Two behaviours appear. Every second- and third-order transfer decays monotonically from
the first step, ending at 22% to 40% of its initial defect. First-order nearest-neighbour
instead overshoots to 4.7 times its first-step value by step 5 and then decays to the
noise floor, so the response of the least accurate transfer is not monotone in time and
its first-step value understates the excursion that follows.

## 14. A supermesh defect in the quadratic projection

The third-order Galerkin projection requires a correction to the supermesh assembly. This
section records the defect, the test that exposes it, and the fix, because the failure it
produces is silent: it yields a plausible restart with no warning, no solver diagnostic
and no visible mesh problem.

**The test.** A Galerkin L2 projection must reproduce exactly any function already lying
in its own space, and constants lie in the quadratic space. Projecting f = 1 through the
uncorrected quadratic path returns values as large as **11.02**, with 380 of 22 508 nodes
wrong by more than 1e-6 and the worst error on the farfield boundary. The linear patch
test f = x fails likewise, with a maximum error of 658. The same tests through the linear
path are satisfied to 5.9e-12.

**The cause.** Coverage of a target cell — its summed intersection measure with the source
cells, divided by its own measure — must be exactly 1 for a conforming tiling, and both
meshes tile the domain exactly. Seven target cells in the outermost ring, at radius 93.2
to 97.2, instead reach a coverage of up to **2.227**.

The clipper can return, for particular configurations, a vertex list that is a valid
*signed* traversal but not a simple polygon: it carries a duplicated vertex and a vertex
lying outside the intersection, so the traversal doubles back on itself. Its shoelace area
is exact, because the spurious lobes carry opposite sign and cancel. The quadratic
assembler, however, fan-triangulated that list and took the **absolute value** of each
sub-triangle area, which adds the lobes instead of cancelling them. For the worst cell the
true intersection area is 6.195763, the shoelace form recovers exactly that, and the
absolute-value fan returns 26.943038; combined with the cell's other contribution the
coverage becomes (10.716094 + 26.943038)/16.911857 = 2.2268, matching the measured value.

The same inflated measure scales the quadrature weights, so the right-hand side is
over-assembled by that factor while the mass matrix, built from the true cell measure, is
not. The resulting solve carries degrees of freedom as large as 27.6 for a density field
bounded in [0.987, 1.008], and the error spreads inward from the outer ring.

Everything else in the path is sound. The quadratic source reconstruction is well behaved,
its edge-midpoint correction being 1% to 3% of the field jump across the edge in every
band with a far-field maximum of 1.1e-7, because the Hessian decays faster than the
squared edge length grows. The mass system solves to a relative residual of 5.4e-14 in 30
iterations and agrees with a direct sparse solve to 3.2e-12. No cell is under-covered, and
the cells carrying the bad degrees of freedom have an edge-ratio of 1.22 against a mesh
median of 1.19.

**Why the consequences are global.** The Stage-2 limiter has unbounded spatial reach. It
clips each node into the range of the source field over the cells that fed it, and because
clipping changes the total it then restores the lumped nodal integral by redistributing the
clip defect over every node that still has room, in proportion to that room. The
redistribution is global and the integral is volume-weighted, so the defect settles
preferentially into the largest control volumes.

Applied to a corrupted projection this converts a failure in seven cells into a bias over
the whole field. The two configurations separate the roles cleanly: without the limiter the
uncorrected projection matches the second-order variant in every band except the outermost,
where its error is 8.24e+01; with the limiter that band reads 6.18e-3 while every other band
carries roughly an order of magnitude more error than the second-order variant. The limiter
is doing what it is designed to do — bounding the values and conserving the lumped integral
— but its repair spreads a local failure rather than isolating it, and that behaviour is
independent of the defect below.

**The correction.** The fan must take the signed sub-triangle area. The linear assembler
already measures its intersection polygons through a signed area and moment formula and is
therefore immune to the same degenerate output, so this brings the quadratic assembler into
line with the convention already in use.

With it, coverage is exactly 1 on every cell and the constant and linear patch tests are
satisfied to 3.9e-11 and 1.1e-10. The results reported throughout Part II are those of the
corrected assembler: a defect norm of 7.126e-4, the smallest of the twelve transfers, a
density integral defect of −7.853e-4 in line with every other method, 14 limiter
activations against 11 for the second-order variant, a maximum surface |ΔC_p| of 0.30, and
a restart response that decays from the first step.

The three-dimensional assembler must **not** receive the same change. There the pieces are
disjoint sub-tetrahedra produced by plane clipping rather than a fan over a single
traversal, so their measures add and the absolute value is correct; taking signed volumes
would subtract legitimate pieces.

**Guarding it.** Two checks are added, neither of which fixes the defect but either of
which would have exposed it. A coverage test that fails outright when any cell exceeds 1,
placed outside the adaptive refinement loop — that loop answers under-coverage by widening
the candidate search, which for an over-covered cell only adds further spurious
contributions. Its threshold is 1e-6 rather than the 1e-9 used for under-coverage:
summing many intersection measures leaves a round-off excess of order 1e-9 on a well-formed
supermesh, and a uniform 100 000-cell case peaks at 1 + 5.3e-9, whereas a genuine
double-count is a first-order error. And a constant-field patch test, which is
mesh-independent, needs no reference data, and detects this entire class of assembly error.

The BDF2 relations held on the uncorrected transfer as precisely as on any other, so the
defect was confined to the transfer and never touched the restart algebra applied to it.

---

<!-- RANS-DATA-BEGIN -->

# Part III — RANS extension: data for the results section

Same experiment as Part I on a wall-resolved hybrid grid with the Spalart–Allmaras (SA) and Menter SST turbulence models. Two additions against the inviscid study: the transported turbulence state is part of what is transferred, and the drag is reported split into its pressure and friction parts. This part collects the measured quantities; it does not interpret them.

## III.1 Setup

**Grid.** `Meshes/Coarse_Sharp.su2`: 23271 nodes, 9686 quadrilaterals in the wall layer and wake block and 26851 triangles outside them; 256 airfoil nodes (256 surface elements), farfield at 100 chords. First off-wall spacing 5.00e-06 to 5.25e-06; the quadrilateral layer extends to 0.110 chords from the wall. The trailing edge is closed and sharp (one node at x = 1). The same grid, undeformed, is used by both models and by the references.

**SA solver settings** (`Base_Interp_RANS_SA/Unsteady_SecondOrder.cfg`): SOLVER = RANS, KIND_TURB_MODEL = SA, TIME_MARCHING = DUAL_TIME_STEPPING-2ND_ORDER, TIME_STEP = 1e-3, INNER_ITER = 800, MACH_NUMBER = 0.1, AOA = 6.0, REYNOLDS_NUMBER = 6.0E6, REF_DIMENSIONALIZATION = FREESTREAM_VEL_EQ_MACH, NUM_METHOD_GRAD = WEIGHTED_LEAST_SQUARES, CFL_NUMBER = 6.0, LINEAR_SOLVER = FGMRES, LINEAR_SOLVER_PREC = ILU, LINEAR_SOLVER_ITER = 60, CONV_NUM_METHOD_FLOW = HLLC, MUSCL_FLOW = YES, SLOPE_LIMITER_FLOW = VAN_ALBADA_EDGE, CONV_NUM_METHOD_TURB = SCALAR_UPWIND, MUSCL_TURB = NO.

**SST solver settings** (`Base_Interp_RANS_SST/Unsteady_SecondOrder.cfg`): SOLVER = RANS, KIND_TURB_MODEL = SST, TIME_MARCHING = DUAL_TIME_STEPPING-2ND_ORDER, TIME_STEP = 1e-3, INNER_ITER = 800, MACH_NUMBER = 0.1, AOA = 6.0, REYNOLDS_NUMBER = 6.0E6, REF_DIMENSIONALIZATION = FREESTREAM_VEL_EQ_MACH, NUM_METHOD_GRAD = WEIGHTED_LEAST_SQUARES, CFL_NUMBER = 6.0, LINEAR_SOLVER = FGMRES, LINEAR_SOLVER_PREC = ILU, LINEAR_SOLVER_ITER = 60, CONV_NUM_METHOD_FLOW = HLLC, MUSCL_FLOW = YES, SLOPE_LIMITER_FLOW = VAN_ALBADA_EDGE, CONV_NUM_METHOD_TURB = SCALAR_UPWIND, MUSCL_TURB = NO.

**Transfer and mesh motion** (both models): InterpMethod = ConsGalerkinProj, InterpOrder = Second, WLSDegree = 1, BoundaryTreatment = False, BoundaryPointValue = False, EnforceConservation = False, ConservationCorrection = True, GalerkinLimit = True, findOneToOneCorrespondence = True, MinWallDistance = 0.01, MaxWallDistance = 20, FreezeBoundary = True, Field = "white", SafeFactor = 4.0, QualityFloor = 0.2, Seed = 12345. The launcher sets the drift band and one-shot bands per case (see the case names). The wall band `MinWallDistance` and every boundary node are held fixed; the wall-layer quadrilaterals therefore never move. The transferred fields are Density, Momentum_x, Momentum_y, Energy and the model's transported turbulence variables (Nu_Tilde for SA; Turb_Kin_Energy and Omega for SST, both transferred as rho-weighted quantities).

**Time units (recovered from the SU2 startup tables of the executed runs).** `REF_DIMENSIONALIZATION = FREESTREAM_VEL_EQ_MACH` with a dimensional free stream (T = 300 K): the reference velocity is the sound speed 347.224 m/s, the reference time L_ref/a = 0.00287999 s, and the configured `TIME_STEP = 0.001` is in seconds ("Unsteady time step provided by the user (s): 0.001"). The internal non-dimensional step is therefore **0.347224** (a = 1, U_inf = 0.1 in code units), as printed in the startup table "| Time Step | 0.001 | 0.00287999 | s | 0.347224 |". Consequences: an acoustic front travels 0.347 chords per step, the predicted force-response lag for a perturbation at radius d is **2.88·d steps** (below one step for every band used), and 400 steps span 138.9 time units = **13.9 chord flow-throughs**. The Euler study (p = rho = T = 1, reference factor 1) advanced 1e-3 per step: 845·d steps and 0.047 chord flow-throughs over 400 steps. The RANS first-order start-up window uses 1e-4 s = 0.0347 non-dimensional.

**Reference runs.** `Unsteady_Ref_RANS_<MODEL>_Long`: 800 steps on the undeformed grid with no transfer, started from the same two restart levels as every case (the converged steady state and the first-order start-up state written by the launcher). Cases are compared to it step by step.

## III.1b Runtime provenance and numerical controls (executed-run evidence)

**Reference quantities and time step, from the SU2 startup tables** (`rans_runtime_units.csv`):

| family | quantity | dimensional | reference | unit | non-dimensional | user time steps (s) | t_ref (s) |
|---|---|---|---|---|---|---|---|
| RANS SA reference | Density | 3.18973 | 3.18973 | kg/m^3 | 1 | 0.001. |  |
| RANS SA reference | Temperature | 300 | 300 | K | 1 | 0.001. |  |
| RANS SA reference | Time Step | 0.001 | 0.00287999 | s | 0.347224 | 0.001. |  |
| RANS SST reference | Density | 3.18973 | 3.18973 | kg/m^3 | 1 | 0.001. |  |
| RANS SST reference | Temperature | 300 | 300 | K | 1 | 0.001. |  |
| RANS SST reference | Time Step | 0.001 | 0.00287999 | s | 0.347224 | 0.001. |  |
| RANS SA drift case (start-up + windows) | Density | 3.18973 | 3.18973 | kg/m^3 | 1 | 0.0001., 0.001. |  |
| RANS SA drift case (start-up + windows) | Temperature | 300 | 300 | K | 1 | 0.0001., 0.001. |  |
| RANS SA drift case (start-up + windows) | Time Step | 0.0001 | 0.00287999 | s | 0.0347224 | 0.0001., 0.001. |  |
| RANS SA drift case (start-up + windows) | Time Step | 0.001 | 0.00287999 | s | 0.347224 | 0.0001., 0.001. |  |
| RANS SST drift case (start-up + windows) | Density | 3.18973 | 3.18973 | kg/m^3 | 1 | 0.0001., 0.001. |  |
| RANS SST drift case (start-up + windows) | Temperature | 300 | 300 | K | 1 | 0.0001., 0.001. |  |
| RANS SST drift case (start-up + windows) | Time Step | 0.0001 | 0.00287999 | s | 0.0347224 | 0.0001., 0.001. |  |
| RANS SST drift case (start-up + windows) | Time Step | 0.001 | 0.00287999 | s | 0.347224 | 0.0001., 0.001. |  |
| Euler reference (Part I) | Density | 1 | 1 | kg/m^3 | 1 | 0.001. |  |
| Euler reference (Part I) | Temperature | 1 | 1 | K | 1 | 0.001. |  |
| Euler reference (Part I) | Time Step | 0.001 | 1 | s | 0.001 | 0.001. |  |

**Code provenance** (`rans_code_provenance.csv`):

| item | value | date |
|---|---|---|
| SU2_CFD binary | 11736224 | 2025-12-10 |
| SU2 source commit | dcee051d8ee7fba29ebe4d5f373c63a5fa198906 | Fri Jun 27 08:21:49 2025 +0200 |
| SU2 working tree dirty files | 0 |  |
| MeshAdaptation commit (transfer / motion code, as built into the case bases) | 35d865ebdea88fc265a62d32bbba2527a187491c | Wed Sep 9 17:57:44 2026 +0200 |
| MeshAdaptation dirty files at build time | 10 |  |
| InterpSolution.py used by the SA cases (sha256) | d0ca621366e755117d61412bc81965ad3a110e0ef25fcd238a087d3c207604af |  |
| InterpSolution.py used by the SST cases (sha256) | d0ca621366e755117d61412bc81965ad3a110e0ef25fcd238a087d3c207604af |  |
| FunScripts/InterpSolution_Fun.py used by the SA cases (sha256) | 5530368696ee491ef6818612cdf669aeb59949697733fe38ea26f759684c3d81 |  |
| FunScripts/InterpSolution_Fun.py used by the SST cases (sha256) | 5530368696ee491ef6818612cdf669aeb59949697733fe38ea26f759684c3d81 |  |
| RandomMeshMovements.py used by the SA cases (sha256) | 74cbb8a4e603dd2f3a59145b887a22c1536172b900ed8c3d0ef69bca0503057e |  |
| RandomMeshMovements.py used by the SST cases (sha256) | 74cbb8a4e603dd2f3a59145b887a22c1536172b900ed8c3d0ef69bca0503057e |  |
| FunScripts/RandomMeshMovements_Fun.py used by the SA cases (sha256) | 25985929d43c869c4b8f5367f6b851d862c5018e7a46f176bfc9ac8f8ee29718 |  |
| FunScripts/RandomMeshMovements_Fun.py used by the SST cases (sha256) | 25985929d43c869c4b8f5367f6b851d862c5018e7a46f176bfc9ac8f8ee29718 |  |
| RunCase_UnsteadyAdaptation.sh used by the SA cases (sha256) | ad21fcdb2e2f67c533731297950c1e626ff402099de9195374e382439e865e45 |  |
| RunCase_UnsteadyAdaptation.sh used by the SST cases (sha256) | ad21fcdb2e2f67c533731297950c1e626ff402099de9195374e382439e865e45 |  |
| FSI.py used by the SA cases (sha256) | 8ed9f85b994789d8a302d58c4bc4f41990e33ba7b431a071ee387f68c1766eb5 |  |
| FSI.py used by the SST cases (sha256) | 8ed9f85b994789d8a302d58c4bc4f41990e33ba7b431a071ee387f68c1766eb5 |  |

**Mesh replay.** The applied displacement of every one of the 80 events is byte-identical across the four seed-1 drift cases (SA/SST × ConsGalerkinProj/NN): 80 of 80 events identical (`rans_mesh_replay.csv`, md5 of the displacement arrays). The motion is generated from the seed and the event index only (`seed = Seed + event`), never from the flow state. The saved `Grid_Velocity` in the volume files is at most 1.8e-32: the grid is stationary between replacements. After each transfer the solver restarts with `DUAL_TIME_STEPPING-2ND_ORDER` from the two transferred levels (log: "second order in time", "Read flow solution from restart_flow_<k-1>/<k>"); the first-order scheme is used only in the start-up window 0→1 with a ten times smaller step. Note for a time-step refinement: the event index enters the seed, so a run with a different step count must replay the archived displacement files by physical event rather than regenerate them.

**Protected region** (`rans_protected_summary.csv`, node sets in `rans_protected_region.csv`):

| item | value |
|---|---|
| nodes total | 23271 |
| moved (first event) | 12602 |
| fixed | 10669 |
| quad-layer nodes | 9981 |
| quad-layer nodes moved | 1782 |
| quad-layer moved: solver wall distance min | 0.0101 |
| quad-layer moved: wall distance max | 0.1097 |
| fixed nodes with solver wall distance < 0.01 | 8590 |
| moved nodes with solver wall distance < 0.01 | 0 |
| moved nodes: solver wall distance min | 0.0100 |
| mask criterion in the motion code | distance to the nearest airfoil NODE > MinWallDistance = 0.01 (not the solver wall distance) |
| boundary nodes frozen | yes (airfoil 370 + farfield 63) |
| max displacement (first event) | 2.389e-01 |

**Mesh quality over the event sequence** (`rans_mesh_quality.csv`, 340 mesh generations): corner-triangle area ratio deformed/undeformed, per-event minimum between 0.209 and 0.574 (quality floor setting 0.2), nodes reduced by the admissibility backoff: 6 at most. The grid is valid at every event.

**Inner-iteration usage per physical step** (`rans_inner_convergence.csv`). SU2 stops the inner loop when log10 rms[Rho] < CONV_RESIDUAL_MINVAL at a check, or at the INNER_ITER cap. The number of inner iterations used, the density residual level reached at the end of each physical step and the linear-solver effort are listed per run.

| run | steps | inner min | inner max | at cap | log10 rms start (min/med/max) | log10 rms end (min/med/max) | linsol iters (min/med/max) | log10 linsol res (min/med/max) | INNER_ITER | CONV_RESIDUAL_MINVAL | CONV_STARTITER |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | 398 | 10 | 799 | 144 | -12.28/-7.39/-6.39 | -12.19/-12.00/-9.12 | 9/10/10 | -8.80/-8.43/-8.00 | 800 | -12 | 10 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | 398 | 10 | 799 | 169 | -12.28/-7.36/-6.33 | -12.19/-12.00/-9.70 | 9/9/10 | -8.82/-8.31/-8.00 | 800 | -12 | 10 |
| SA 0.01-20 NN Drift | 398 | 10 | 799 | 156 | -12.28/-7.35/-6.41 | -12.19/-12.00/-9.21 | 9/10/10 | -8.85/-8.41/-8.00 | 800 | -12 | 10 |
| SA 0.01-20 NN Drift Seed2 | 398 | 10 | 799 | 171 | -12.28/-7.34/-6.24 | -12.19/-12.00/-9.32 | 9/10/10 | -8.84/-8.37/-8.00 | 800 | -12 | 10 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 398 | 10 | 305 | 0 | -12.41/-9.88/-7.18 | -12.99/-12.32/-12.00 | 9/10/12 | -8.80/-8.22/-8.00 | 800 | -12 | 10 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 398 | 10 | 200 | 0 | -12.41/-9.72/-7.06 | -12.73/-12.29/-12.00 | 9/10/12 | -8.78/-8.31/-8.01 | 800 | -12 | 10 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 398 | 10 | 158 | 0 | -12.41/-10.05/-7.36 | -12.75/-12.31/-12.00 | 9/10/11 | -8.78/-8.30/-8.01 | 800 | -12 | 10 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 798 | 10 | 88 | 0 | -12.61/-10.54/-7.62 | -12.74/-12.44/-12.00 | 9/10/12 | -8.78/-8.57/-8.01 | 800 | -12 | 10 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 398 | 10 | 129 | 0 | -12.41/-10.05/-7.50 | -12.85/-12.31/-12.00 | 9/10/12 | -8.78/-8.31/-8.01 | 800 | -12 | 10 |
| SST 0.01-20 ConsGalerkinProj Drift | 398 | 10 | 799 | 333 | -12.06/-7.39/-6.56 | -12.08/-11.10/-10.11 | 8/9/11 | -8.82/-8.30/-8.00 | 800 | -12 | 10 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | 398 | 10 | 799 | 284 | -12.06/-7.37/-6.31 | -12.08/-11.31/-9.46 | 8/9/11 | -8.86/-8.40/-8.00 | 800 | -12 | 10 |
| SST 0.01-20 NN Drift | 398 | 10 | 799 | 323 | -12.06/-7.36/-6.51 | -12.08/-11.24/-9.52 | 9/9/11 | -8.86/-8.36/-8.00 | 800 | -12 | 10 |
| SST 0.01-20 NN Drift Seed2 | 398 | 10 | 799 | 277 | -12.06/-7.35/-6.34 | -12.08/-11.24/-9.34 | 8/9/11 | -8.87/-8.41/-8.00 | 800 | -12 | 10 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 398 | 10 | 270 | 0 | -12.06/-9.88/-7.18 | -12.26/-12.10/-12.00 | 10/10/11 | -8.71/-8.16/-8.01 | 800 | -12 | 10 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 398 | 10 | 218 | 0 | -12.06/-9.74/-7.07 | -12.24/-12.10/-12.00 | 10/11/12 | -8.75/-8.25/-8.00 | 800 | -12 | 10 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 398 | 10 | 185 | 0 | -12.06/-10.05/-7.37 | -12.22/-12.10/-12.00 | 10/11/11 | -8.76/-8.24/-8.03 | 800 | -12 | 10 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 798 | 10 | 91 | 0 | -12.06/-10.65/-7.64 | -12.34/-12.15/-12.00 | 9/10/12 | -8.73/-8.52/-8.00 | 800 | -12 | 10 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 398 | 10 | 151 | 0 | -12.06/-10.07/-7.51 | -12.24/-12.10/-12.00 | 10/10/12 | -8.78/-8.31/-8.01 | 800 | -12 | 10 |
| Unsteady_Ref_RANS_SA_Long | 798 | 10 | 10 | 0 | -12.78/-12.50/-12.11 | -12.81/-12.54/-12.12 | 9/10/10 | -8.78/-8.43/-8.01 | 800 | -12 | 10 |
| Unsteady_Ref_RANS_SST_Long | 798 | 10 | 10 | 0 | -12.08/-12.03/-11.34 | -12.27/-12.17/-12.06 | 8/9/11 | -8.72/-8.24/-8.00 | 800 | -12 | 10 |

## III.2 Steady solutions and references

Converged steady state (`Steady_RANS_<MODEL>/forces_breakdown.dat`, limiter on, 30 000 iterations, residual of density below 1e-8):

| model | CL | CD | CD pressure | CD friction | CMz |
|---|---|---|---|---|---|
| SA | 0.676257 | 0.010229 | 0.003559 | 0.006671 | 0.000671 |
| SST | 0.675605 | 0.009858 | 0.003469 | 0.006389 | 0.000382 |

Reference (no transfer) force histories:

| model | steps | CL first | CL step 399 | CL last | CD first | CD last | CL std | CD std |
|---|---|---|---|---|---|---|---|---|
| SA | 2..799 | 0.676257 | 0.676257 | 0.676257 | 0.010229 | 0.010229 | 1.317e-09 | 9.180e-10 |
| SST | 2..799 | 0.675605 | 0.675603 | 0.675603 | 0.009858 | 0.009858 | 3.749e-07 | 2.391e-08 |

## III.2b Reference on a moved grid

The final grid of the SA ConsGalerkinProj drift case after its 40 events (`Meshes/Coarse_Sharp_moved_00399.su2`; the mesh sequence is seed-determined and identical for both models; 12602 of 23271 nodes move at an event, wall layer and boundary fixed) is run without any transfer: the steady RANS converged on it from the original-grid steady state, then 400 unsteady steps from that state. Differences to the same runs on the original grid are the effect of the grid alone. The drift cases compute their first window, steps 10..19, on the first moved grid from the transferred state; their deviation in that window is listed for comparison (`rans_deformed_reference.csv`, `RANS_<MODEL>_DeformedReference.png`).

| model | what | n | dCL | dCD | dCDp | dCDv | dCMz |
|---|---|---|---|---|---|---|---|
| SA | steady | 19999 | 6.912e-04 | 1.054e-04 | 1.040e-04 | 1.000e-06 | 1.834e-04 |
| SA | unsteady 100..399 | 300 | 7.002e-04 (trend -9.9e-07 ± 5.2e-07) | 1.059e-04 (trend -1.5e-07 ± 1.1e-07) |  |  | 1.848e-04 (trend -2.5e-07 ± 2.8e-07) |
| SA | drift case first window 10..19: RANS_SA_WD_0.01-20_ConsGalerkinProj_Drift_Second | 10 | 2.788e-05 (min 1.4e-05, max 5.5e-05) | 3.778e-05 (min 3.1e-05, max 4.8e-05) |  |  | 2.818e-06 (min -1.8e-05, max 2.4e-05) |
| SA | drift case first window 10..19: RANS_SA_WD_0.01-20_NN_Drift_Second | 10 | 2.949e-05 (min 1.6e-05, max 5.2e-05) | 3.791e-05 (min 3.2e-05, max 4.9e-05) |  |  | 3.267e-06 (min -1.7e-05, max 2.4e-05) |
| SA | steady, surface integration | 0 | 6.910e-04 | 1.053e-04 | 1.043e-04 | 1.087e-06 |  |
| SST | steady | 19999 | 6.249e-04 | 8.671e-05 | 8.800e-05 | -1.000e-06 | 1.852e-04 |
| SST | unsteady 100..399 | 300 | 6.169e-04 (trend -2.1e-08 ± 1.5e-08) | 8.574e-05 (trend 3.8e-09 ± 1.2e-08) |  |  | 1.810e-04 (trend 7.2e-09 ± 3.2e-08) |
| SST | drift case first window 10..19: RANS_SST_WD_0.01-20_ConsGalerkinProj_Drift_Second | 10 | -3.706e-05 (min -6.4e-05, max 1.9e-05) | 3.548e-05 (min 2.7e-05, max 4.6e-05) |  |  | -2.633e-05 (min -5.6e-05, max 1.7e-07) |
| SST | drift case first window 10..19: RANS_SST_WD_0.01-20_NN_Drift_Second | 10 | -3.500e-05 (min -6.3e-05, max 1.6e-05) | 3.569e-05 (min 2.8e-05, max 4.4e-05) |  |  | -2.601e-05 (min -5.4e-05, max -1.7e-06) |
| SST | steady, surface integration | 0 | 6.249e-04 | 8.671e-05 | 8.744e-05 | -7.369e-07 |  |

Surface difference between the two steady states, node to node on the 256 airfoil nodes:

| model | RMS dCp | max dCp | x of max | RMS dCf_x | max dCf_x | x of max |
|---|---|---|---|---|---|---|
| SA | 3.107e-03 | 3.595e-02 | 1.0000 | 1.042e-05 | 6.326e-05 | 0.2011 |
| SST | 2.230e-03 | 1.160e-02 | 0.7473 | 1.083e-05 | 6.740e-05 | 0.2011 |

## III.3 Force drift under repeated transfer (transfer every 10 steps, 400 steps)

Deviation from the reference, OLS trend per 100 steps with the standard error of the slope. CL, CD, CMz from the history files (full precision); CDp, CDv from `forces_breakdown` (6 decimals; "rounding" marks a signal below that precision). The OLS standard error assumes independent samples; the residuals have a lag-1 autocorrelation of 0.8 to 0.95, so an AR(1)-corrected error se·sqrt((1+r1)/(1−r1)) and the effective sample size N·(1−r1)/(1+r1) are listed too. "resolved" = |trend| > 2 standard errors, once with the OLS error and once with the AR(1) error. All slopes are descriptive finite-window fits over one 400-step sequence; events within a sequence are correlated observations, not independent experiments.

| case | coeff | ref mean | offset | scatter | max |dev| | trend /100it | ± OLS | resolved OLS | lag-1 r | ± AR(1) | N_eff | resolved AR(1) | % /100it |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | CL | 6.762569e-01 | 7.232e-05 | 2.098e-04 | 8.266e-04 | -6.742e-05 | 8.53e-06 | yes | 0.939 | 4.82e-05 | 12 | no | -0.0100 |
| SA 0.01-20 ConsGalerkinProj Drift | CD | 1.022935e-02 | 5.923e-05 | 2.544e-05 | 1.269e-04 | -1.252e-06 | 1.11e-06 | no | 0.930 | 5.84e-06 | 14 | no | -0.0122 |
| SA 0.01-20 ConsGalerkinProj Drift | CDp | 3.559000e-03 | 5.684e-05 | 2.556e-05 | 1.250e-04 | -1.345e-06 | 1.12e-06 | no | 0.931 | 5.88e-06 | 14 | no | -0.0378 |
| SA 0.01-20 ConsGalerkinProj Drift | CDv | 6.671000e-03 | 1.754e-06 | 1.418e-06 | 6.000e-06 | 6.331e-08 | 6.19e-08 | no | 0.970 | 5.01e-07 | 6 | no | 0.0009 |
| SA 0.01-20 ConsGalerkinProj Drift | CMz | 6.708957e-04 | 2.747e-05 | 9.099e-05 | 3.358e-04 | -1.892e-05 | 3.86e-06 | yes | 0.919 | 1.89e-05 | 17 | no | -2.8196 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | CL | 6.762569e-01 | -3.523e-05 | 2.148e-04 | 6.303e-04 | -2.051e-05 | 9.34e-06 | yes | 0.940 | 5.30e-05 | 12 | no | -0.0030 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | CD | 1.022935e-02 | 5.580e-05 | 2.397e-05 | 1.200e-04 | 3.983e-06 | 1.03e-06 | yes | 0.919 | 5.01e-06 | 17 | no | 0.0389 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | CDp | 3.559000e-03 | 5.430e-05 | 2.381e-05 | 1.170e-04 | 4.274e-06 | 1.02e-06 | yes | 0.917 | 4.90e-06 | 17 | no | 0.1201 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | CDv | 6.671000e-03 | 8.317e-07 | 2.563e-06 | 1.300e-05 | -3.362e-07 | 1.11e-07 | yes | 0.986 | 1.34e-06 | 3 | no | -0.0050 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | CMz | 6.708957e-04 | 9.099e-06 | 9.410e-05 | 2.779e-04 | -2.344e-06 | 4.11e-06 | no | 0.919 | 2.00e-05 | 17 | no | -0.3494 |
| SA 0.01-20 NN Drift | CL | 6.762569e-01 | 5.929e-05 | 2.111e-04 | 7.782e-04 | -6.585e-05 | 8.62e-06 | yes | 0.943 | 5.02e-05 | 12 | no | -0.0097 |
| SA 0.01-20 NN Drift | CD | 1.022935e-02 | 5.819e-05 | 2.587e-05 | 1.241e-04 | -6.967e-07 | 1.13e-06 | no | 0.935 | 6.16e-06 | 13 | no | -0.0068 |
| SA 0.01-20 NN Drift | CDp | 3.559000e-03 | 5.591e-05 | 2.591e-05 | 1.220e-04 | -7.605e-07 | 1.13e-06 | no | 0.935 | 6.18e-06 | 13 | no | -0.0214 |
| SA 0.01-20 NN Drift | CDv | 6.671000e-03 | 1.646e-06 | 1.424e-06 | 6.000e-06 | 2.040e-08 | 6.23e-08 | no | 0.969 | 4.96e-07 | 6 | no | 0.0003 |
| SA 0.01-20 NN Drift | CMz | 6.708957e-04 | 2.489e-05 | 9.139e-05 | 3.258e-04 | -1.807e-05 | 3.89e-06 | yes | 0.923 | 1.94e-05 | 16 | no | -2.6933 |
| SA 0.01-20 NN Drift Seed2 | CL | 6.762569e-01 | -2.451e-05 | 2.117e-04 | 6.536e-04 | -1.511e-05 | 9.23e-06 | no | 0.941 | 5.29e-05 | 12 | no | -0.0022 |
| SA 0.01-20 NN Drift Seed2 | CD | 1.022935e-02 | 5.621e-05 | 2.331e-05 | 1.167e-04 | 3.790e-06 | 1.00e-06 | yes | 0.918 | 4.84e-06 | 17 | no | 0.0371 |
| SA 0.01-20 NN Drift Seed2 | CDp | 3.559000e-03 | 5.461e-05 | 2.321e-05 | 1.160e-04 | 4.190e-06 | 9.93e-07 | yes | 0.916 | 4.74e-06 | 17 | no | 0.1177 |
| SA 0.01-20 NN Drift Seed2 | CDv | 6.671000e-03 | 9.271e-07 | 2.869e-06 | 1.400e-05 | -4.292e-07 | 1.24e-07 | yes | 0.988 | 1.59e-06 | 2 | no | -0.0064 |
| SA 0.01-20 NN Drift Seed2 | CMz | 6.708957e-04 | 1.183e-05 | 9.401e-05 | 2.734e-04 | 4.530e-08 | 4.11e-06 | no | 0.920 | 2.02e-05 | 17 | no | 0.0068 |
| SST 0.01-20 ConsGalerkinProj Drift | CL | 6.756033e-01 | -9.986e-05 | 2.423e-04 | 8.247e-04 | -7.922e-05 | 9.82e-06 | yes | 0.935 | 5.37e-05 | 13 | no | -0.0117 |
| SST 0.01-20 ConsGalerkinProj Drift | CD | 9.858431e-03 | 6.000e-05 | 2.992e-05 | 1.369e-04 | 1.864e-06 | 1.31e-06 | no | 0.929 | 6.80e-06 | 15 | no | 0.0189 |
| SST 0.01-20 ConsGalerkinProj Drift | CDp | 3.469000e-03 | 5.997e-05 | 2.987e-05 | 1.360e-04 | 1.809e-06 | 1.30e-06 | no | 0.929 | 6.79e-06 | 15 | no | 0.0521 |
| SST 0.01-20 ConsGalerkinProj Drift | CDv | 6.389000e-03 | 4.221e-07 | 1.231e-06 | 3.000e-06 | 1.247e-07 | 5.35e-08 | yes (rounding) | 0.961 | 3.81e-07 | 8 | no | 0.0020 |
| SST 0.01-20 ConsGalerkinProj Drift | CMz | 3.819379e-04 | -4.037e-06 | 1.050e-04 | 3.285e-04 | -1.779e-05 | 4.50e-06 | yes | 0.918 | 2.18e-05 | 17 | no | -4.6587 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | CL | 6.756033e-01 | -2.035e-04 | 2.748e-04 | 1.010e-03 | -7.879e-05 | 1.13e-05 | yes | 0.939 | 6.39e-05 | 13 | no | -0.0117 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | CD | 9.858431e-03 | 5.845e-05 | 3.072e-05 | 1.288e-04 | 3.075e-06 | 1.33e-06 | yes | 0.927 | 6.88e-06 | 15 | no | 0.0312 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | CDp | 3.469000e-03 | 5.858e-05 | 3.086e-05 | 1.290e-04 | 2.916e-06 | 1.34e-06 | yes | 0.927 | 6.92e-06 | 15 | no | 0.0841 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | CDv | 6.389000e-03 | 3.040e-07 | 2.697e-06 | 1.400e-05 | 1.862e-07 | 1.18e-07 | no | 0.979 | 1.15e-06 | 4 | no | 0.0029 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | CMz | 3.819379e-04 | -2.311e-05 | 1.153e-04 | 3.458e-04 | -1.512e-05 | 4.99e-06 | yes | 0.921 | 2.47e-05 | 16 | no | -3.9598 |
| SST 0.01-20 NN Drift | CL | 6.756033e-01 | -9.703e-05 | 2.456e-04 | 8.222e-04 | -7.924e-05 | 9.98e-06 | yes | 0.935 | 5.45e-05 | 13 | no | -0.0117 |
| SST 0.01-20 NN Drift | CD | 9.858431e-03 | 6.069e-05 | 3.041e-05 | 1.375e-04 | 1.730e-06 | 1.33e-06 | no | 0.928 | 6.85e-06 | 15 | no | 0.0176 |
| SST 0.01-20 NN Drift | CDp | 3.469000e-03 | 6.054e-05 | 3.042e-05 | 1.370e-04 | 1.667e-06 | 1.33e-06 | no | 0.928 | 6.85e-06 | 15 | no | 0.0480 |
| SST 0.01-20 NN Drift | CDv | 6.389000e-03 | 5.628e-07 | 1.311e-06 | 4.000e-06 | 5.680e-08 | 5.73e-08 | no (rounding) | 0.962 | 4.11e-07 | 8 | no | 0.0009 |
| SST 0.01-20 NN Drift | CMz | 3.819379e-04 | -1.994e-06 | 1.067e-04 | 3.230e-04 | -1.711e-05 | 4.59e-06 | yes | 0.918 | 2.22e-05 | 17 | no | -4.4789 |
| SST 0.01-20 NN Drift Seed2 | CL | 6.756033e-01 | -1.963e-04 | 2.530e-04 | 9.459e-04 | -5.603e-05 | 1.07e-05 | yes | 0.936 | 5.89e-05 | 13 | no | -0.0083 |
| SST 0.01-20 NN Drift Seed2 | CD | 9.858431e-03 | 5.937e-05 | 2.990e-05 | 1.311e-04 | 4.743e-06 | 1.29e-06 | yes | 0.930 | 6.77e-06 | 14 | no | 0.0481 |
| SST 0.01-20 NN Drift Seed2 | CDp | 3.469000e-03 | 6.007e-05 | 2.976e-05 | 1.300e-04 | 5.068e-06 | 1.28e-06 | yes | 0.929 | 6.64e-06 | 15 | no | 0.1461 |
| SST 0.01-20 NN Drift Seed2 | CDv | 6.389000e-03 | -3.116e-07 | 2.026e-06 | 8.000e-06 | -2.996e-07 | 8.73e-08 | yes | 0.979 | 8.45e-07 | 4 | no | -0.0047 |
| SST 0.01-20 NN Drift Seed2 | CMz | 3.819379e-04 | -1.881e-05 | 1.086e-04 | 3.446e-04 | -7.450e-06 | 4.74e-06 | no | 0.919 | 2.30e-05 | 17 | no | -1.9505 |

Drag split at full precision (surface integration of Cp and Cf, `rans_surface_forces.csv`), same OLS trend per 100 steps:

| case | steps | dCL mean | dCL max | dCL trend | dCDp mean | dCDp max | dCDp trend | dCDv mean | dCDv max | dCDv trend |
|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | 478 | 7.515e-05 | 8.266e-04 | -6.758e-05 ± 8.0e-06 * | 5.794e-05 | 1.253e-04 | -1.170e-06 ± 1.0e-06 | 1.930e-06 | 6.446e-06 | 6.570e-08 ± 5.6e-08 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | 478 | -3.033e-05 | 6.302e-04 | -2.054e-05 ± 8.8e-06 * | 5.552e-05 | 1.173e-04 | 4.504e-06 ± 9.6e-07 * | 1.028e-06 | 1.287e-05 | -3.048e-07 ± 9.9e-08 * |
| SA 0.01-20 NN Drift | 478 | 6.201e-05 | 7.783e-04 | -6.574e-05 ± 8.0e-06 * | 5.680e-05 | 1.223e-04 | -4.904e-07 ± 1.0e-06 | 1.840e-06 | 6.069e-06 | 3.833e-08 ± 5.5e-08 |
| SA 0.01-20 NN Drift Seed2 | 478 | -2.077e-05 | 6.536e-04 | -1.538e-05 ± 8.7e-06 | 5.552e-05 | 1.161e-04 | 4.363e-06 ± 9.3e-07 * | 1.141e-06 | 1.453e-05 | -4.023e-07 ± 1.1e-07 * |
| SST 0.01-20 ConsGalerkinProj Drift | 478 | -9.673e-05 | 8.246e-04 | -7.858e-05 ± 9.0e-06 * | 6.044e-05 | 1.361e-04 | 1.926e-06 ± 1.2e-06 | 2.835e-07 | 2.968e-06 | 9.766e-08 ± 4.7e-08 * |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | 478 | -2.002e-04 | 1.010e-03 | -8.123e-05 ± 1.1e-05 * | 5.912e-05 | 1.285e-04 | 2.513e-06 ± 1.2e-06 * | 8.818e-08 | 1.399e-05 | 1.727e-07 ± 1.1e-07 |
| SST 0.01-20 NN Drift | 478 | -9.476e-05 | 8.222e-04 | -7.798e-05 ± 9.2e-06 * | 6.076e-05 | 1.367e-04 | 1.928e-06 ± 1.2e-06 | 3.982e-07 | 3.330e-06 | 6.495e-08 ± 4.9e-08 |
| SST 0.01-20 NN Drift Seed2 | 478 | -1.925e-04 | 9.458e-04 | -5.850e-05 ± 9.9e-06 * | 6.069e-05 | 1.296e-04 | 4.688e-06 ± 1.2e-06 * | -4.902e-07 | 7.709e-06 | -2.929e-07 ± 7.9e-08 * |

(* = resolved at 2 standard errors)

## III.4 Response to a single transfer

One transfer at the event step, then as many steps again. Onset lag = first step after the event at which |deviation| reaches the given fraction of that case's own peak. "pre-floor" is the largest deviation before the event (bit-identity check); "leak" the deviation at the event step. CENSORED = the peak lies in the last 10% of the window.

**CL**

| case | d | event | window | pre-floor | leak | peak | peak lag | onset 1% | onset 5% | onset 10% | onset 25% | onset 50% | predicted lag | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 2.31e-05 | 9.849e-05 | 199 | 0 | 0 | 0 | 1 | 1 | 0.14 | CENSORED |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 4.33e-05 | 1.276e-04 | 199 | 0 | 0 | 0 | 0 | 1 | 0.29 | CENSORED |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 2.04e-05 | 6.653e-05 | 199 | 0 | 0 | 0 | 0 | 1 | 0.29 | CENSORED |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 1.06e-05 | 2.581e-05 | 399 | 0 | 0 | 0 | 0 | 1 | 0.58 | CENSORED |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 3.14e-06 | 5.005e-05 | 199 | 0 | 0 | 1 | 2 | 4 | 0.58 | CENSORED |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 2.38e-05 | 1.158e-04 | 199 | 0 | 0 | 0 | 1 | 1 | 0.14 | CENSORED |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 4.13e-05 | 1.229e-04 | 199 | 0 | 0 | 0 | 0 | 1 | 0.29 | CENSORED |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 2.05e-05 | 6.918e-05 | 199 | 0 | 0 | 0 | 0 | 2 | 0.29 | CENSORED |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 9.95e-06 | 2.202e-05 | 8 | 0 | 0 | 0 | 0 | 1 | 0.58 |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 3.29e-06 | 4.748e-05 | 199 | 0 | 0 | 1 | 2 | 4 | 0.58 | CENSORED |

**CDp**

| case | d | event | window | pre-floor | leak | peak | peak lag | onset 1% | onset 5% | onset 10% | onset 25% | onset 50% | predicted lag | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 1.00e-06 | 9.000e-06 | 3 | 0 | 0 | 0 | 1 | 1 | 0.14 | rounding |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 2.00e-06 | 1.100e-05 | 2 | 0 | 0 | 0 | 1 | 1 | 0.29 | rounding |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 2.00e-06 | 5.000e-06 | 2 | 0 | 0 | 0 | 0 | 1 | 0.29 | rounding |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 2.00e-06 | 3.000e-06 | 1 | 0 | 0 | 0 | 0 | 0 | 0.58 | rounding |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 2.00e-06 | 4.000e-06 | 4 | 0 | 0 | 0 | 0 | 0 | 0.58 | rounding |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 0.00e+00 | 1.000e-05 | 3 | 1 | 1 | 1 | 1 | 1 | 0.14 | rounding |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 3.00e-06 | 1.100e-05 | 2 | 0 | 0 | 0 | 0 | 1 | 0.29 | rounding |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 3.00e-06 | 6.000e-06 | 2 | 0 | 0 | 0 | 0 | 0 | 0.29 | rounding |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 1.00e-06 | 2.000e-06 | 1 | 0 | 0 | 0 | 0 | 1 | 0.58 | rounding |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 1.00e-06 | 3.000e-06 | 5 | 0 | 0 | 0 | 0 | 1 | 0.58 | rounding |

**CDv**

| case | d | event | window | pre-floor | leak | peak | peak lag | onset 1% | onset 5% | onset 10% | onset 25% | onset 50% | predicted lag | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 0.14 | rounding |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 0.29 | rounding |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 0.29 | rounding |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 0.00e+00 | 1.000e-06 | 169 | 169 | 169 | 169 | 169 | 169 | 0.58 | rounding |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 0.58 | rounding |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 0.14 | rounding |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 0.29 | rounding |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 0.29 | rounding |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 0.58 | rounding |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 0.58 | rounding |

**CD**

| case | d | event | window | pre-floor | leak | peak | peak lag | onset 1% | onset 5% | onset 10% | onset 25% | onset 50% | predicted lag | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 4.85e-07 | 9.276e-06 | 4 | 0 | 0 | 1 | 1 | 1 | 0.14 |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 2.65e-06 | 1.172e-05 | 2 | 0 | 0 | 0 | 1 | 1 | 0.29 |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 2.19e-06 | 5.536e-06 | 2 | 0 | 0 | 0 | 0 | 1 | 0.29 |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 1.85e-06 | 2.631e-06 | 1 | 0 | 0 | 0 | 0 | 0 | 0.58 |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 1.28e-06 | 3.855e-06 | 6 | 0 | 0 | 0 | 0 | 1 | 0.58 |  |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 1.73e-07 | 9.955e-06 | 4 | 0 | 1 | 1 | 1 | 1 | 0.14 |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 2.50e-06 | 1.125e-05 | 2 | 0 | 0 | 0 | 1 | 1 | 0.29 |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 2.32e-06 | 5.614e-06 | 2 | 0 | 0 | 0 | 0 | 1 | 0.29 |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 1.68e-06 | 2.288e-06 | 1 | 0 | 0 | 0 | 0 | 0 | 0.58 |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 1.27e-06 | 3.293e-06 | 6 | 0 | 0 | 0 | 0 | 1 | 0.58 |  |

**CMz**

| case | d | event | window | pre-floor | leak | peak | peak lag | onset 1% | onset 5% | onset 10% | onset 25% | onset 50% | predicted lag | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 1.01e-05 | 2.612e-05 | 3 | 0 | 0 | 0 | 0 | 1 | 0.14 |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 1.38e-05 | 3.078e-05 | 2 | 0 | 0 | 0 | 0 | 1 | 0.29 |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 5.31e-06 | 1.253e-05 | 2 | 0 | 0 | 0 | 0 | 1 | 0.29 |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 1.64e-06 | 3.203e-06 | 2 | 0 | 0 | 0 | 0 | 0 | 0.58 |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 1.30e-06 | 6.983e-06 | 4 | 0 | 0 | 0 | 2 | 3 | 0.58 |  |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 1.05e-05 | 2.755e-05 | 3 | 0 | 0 | 0 | 0 | 1 | 0.14 |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 1.26e-05 | 2.940e-05 | 2 | 0 | 0 | 0 | 0 | 1 | 0.29 |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 5.39e-06 | 1.280e-05 | 2 | 0 | 0 | 0 | 0 | 1 | 0.29 |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 1.35e-06 | 2.286e-06 | 1 | 0 | 0 | 0 | 0 | 0 | 0.58 |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 1.16e-06 | 6.430e-06 | 4 | 0 | 0 | 0 | 2 | 3 | 0.58 |  |

Lag scaling fit lag = C·d + b over the event-at-200 cases (prediction C = 2.88):

| model | coeff | threshold | cases | C | b | R² | C / predicted |
|---|---|---|---|---|---|---|---|
| SA | CL | 0.01 | 3 | 0.0 | 0.00 | nan | 0.000 |
| SA | CL | 0.05 | 3 | 0.0 | 0.00 | nan | 0.000 |
| SA | CL | 0.1 | 3 | 7.1 | -0.50 | 0.8929 | 2.480 |
| SA | CL | 0.25 | 3 | 8.6 | -0.00 | 0.4286 | 2.976 |
| SA | CL | 0.5 | 3 | 21.4 | -0.50 | 0.8929 | 7.441 |
| SA | CDp | 0.01 | 3 | 0.0 | 0.00 | nan | 0.000 |
| SA | CDp | 0.05 | 3 | 0.0 | 0.00 | nan | 0.000 |
| SA | CDp | 0.1 | 3 | 0.0 | 0.00 | nan | 0.000 |
| SA | CDp | 0.25 | 3 | -5.7 | 1.00 | 0.5714 | -1.984 |
| SA | CDp | 0.5 | 3 | -7.1 | 1.50 | 0.8929 | -2.480 |
| SST | CL | 0.01 | 3 | 0.0 | 0.00 | nan | 0.000 |
| SST | CL | 0.05 | 3 | 0.0 | 0.00 | nan | 0.000 |
| SST | CL | 0.1 | 3 | 7.1 | -0.50 | 0.8929 | 2.480 |
| SST | CL | 0.25 | 3 | 8.6 | -0.00 | 0.4286 | 2.976 |
| SST | CL | 0.5 | 3 | 20.0 | -0.00 | 1.0000 | 6.944 |
| SST | CDp | 0.01 | 3 | -5.7 | 1.00 | 0.5714 | -1.984 |
| SST | CDp | 0.05 | 3 | -5.7 | 1.00 | 0.5714 | -1.984 |
| SST | CDp | 0.1 | 3 | -5.7 | 1.00 | 0.5714 | -1.984 |
| SST | CDp | 0.25 | 3 | -5.7 | 1.00 | 0.5714 | -1.984 |
| SST | CDp | 0.5 | 3 | 1.4 | 0.50 | 0.0357 | 0.496 |

**SA amplitude control** (d = 0.10, SafeFactor 2 against 4): CL peak ratio 1.92; onset shift in steps by threshold 1%: +0, 5%: +0, 10%: +0, 25%: +0, 50%: +0.

**SST amplitude control** (d = 0.10, SafeFactor 2 against 4): CL peak ratio 1.78; onset shift in steps by threshold 1%: +0, 5%: +0, 10%: +0, 25%: +0, 50%: -1.

Peak deviations after the event from the surface integration (full-precision split):

| case | dCL pre-floor | dCL peak | at lag | dCDp pre-floor | dCDp peak | at lag | dCDv pre-floor | dCDv peak | at lag |
|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.000e+00 | 9.849e-05 | 199 | 0.000e+00 | 9.282e-06 | 4 | 0.000e+00 | -2.941e-07 | 199 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.000e+00 | 1.276e-04 | 199 | 0.000e+00 | 1.172e-05 | 2 | 0.000e+00 | 2.247e-07 | 199 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.000e+00 | 6.653e-05 | 199 | 0.000e+00 | 5.534e-06 | 2 | 0.000e+00 | 1.522e-07 | 199 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.000e+00 | -2.580e-05 | 399 | 0.000e+00 | -2.631e-06 | 1 | 0.000e+00 | -7.376e-07 | 399 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.000e+00 | -5.005e-05 | 199 | 0.000e+00 | -3.855e-06 | 6 | 0.000e+00 | -2.996e-07 | 199 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.000e+00 | 1.158e-04 | 199 | 0.000e+00 | 9.959e-06 | 4 | 0.000e+00 | -2.083e-07 | 199 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.000e+00 | 1.229e-04 | 199 | 0.000e+00 | 1.125e-05 | 2 | 0.000e+00 | 1.692e-07 | 199 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.000e+00 | 6.918e-05 | 199 | 0.000e+00 | 5.613e-06 | 2 | 0.000e+00 | 1.087e-07 | 199 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.000e+00 | -2.202e-05 | 8 | 0.000e+00 | -2.289e-06 | 1 | 0.000e+00 | -6.887e-07 | 399 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.000e+00 | -4.748e-05 | 199 | 0.000e+00 | -3.293e-06 | 6 | 0.000e+00 | -2.370e-07 | 199 |

## III.4b Pressure probe inside the launch annulus

Shell-RMS of the pressure difference case − reference on shells of wall distance r inside the launch radius d (nodes that never moved), as in Part I §6. Onset lag per shell at several fractions of the shell's own peak; the acoustic arrival at a = 1 with the non-dimensional step 0.347 is 2.88·(d − r) steps (`rans_pressure_probe.csv`, `RANS_<MODEL>_PressureProbe_Pressure.png`).

The volume fields are written in single precision; one ULP of the pressure (p ≈ 0.714) is 5.96e-08. Shells whose peak is below 3 ULP are marked: their onset lags are rounding, not arrivals. The force histories are double precision and are not affected.

| case | d | shell r | nodes | d − r | acoustic lag | peak | onset 2% | onset 5% | onset 10% | onset 25% | onset 50% | onset 75% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.005 | 3139 | 0.045 | 0.13 | 7.745e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.01 | 1118 | 0.040 | 0.12 | 7.166e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.02 | 577 | 0.030 | 0.09 | 7.973e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.04 | 329 | 0.010 | 0.03 | 2.169e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.005 | 3139 | 0.095 | 0.27 | 1.327e-06 | 0 | 0 | 0 | 1 | 4 | 6 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.01 | 1118 | 0.090 | 0.26 | 9.380e-07 | 0 | 0 | 0 | 0 | 3 | 5 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.02 | 577 | 0.080 | 0.23 | 6.621e-07 | 0 | 0 | 0 | 0 | 1 | 4 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.04 | 329 | 0.060 | 0.17 | 6.509e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.06 | 258 | 0.040 | 0.12 | 1.461e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.08 | 201 | 0.020 | 0.06 | 3.452e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.005 | 3139 | 0.095 | 0.27 | 5.738e-07 | 0 | 0 | 0 | 0 | 4 | 6 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.01 | 1118 | 0.090 | 0.26 | 4.176e-07 | 0 | 0 | 0 | 0 | 1 | 5 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.02 | 577 | 0.080 | 0.23 | 3.228e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.04 | 329 | 0.060 | 0.17 | 3.581e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.06 | 258 | 0.040 | 0.12 | 8.044e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.08 | 201 | 0.020 | 0.06 | 1.731e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.005 | 3139 | 0.195 | 0.56 | 4.824e-07 | 0 | 0 | 0 | 4 | 7 | 8 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.01 | 1118 | 0.190 | 0.55 | 4.071e-07 | 0 | 0 | 0 | 4 | 7 | 9 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.02 | 577 | 0.180 | 0.52 | 3.235e-07 | 0 | 0 | 0 | 1 | 7 | 8 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.04 | 329 | 0.160 | 0.46 | 1.733e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 1 | 7 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.06 | 258 | 0.140 | 0.40 | 1.144e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.08 | 201 | 0.120 | 0.35 | 1.011e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.1 | 116 | 0.100 | 0.29 | 1.190e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.13 | 137 | 0.070 | 0.20 | 3.285e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.16 | 88 | 0.040 | 0.12 | 6.304e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.18 | 66 | 0.020 | 0.06 | 4.697e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.005 | 3139 | 0.195 | 0.56 | 5.355e-07 | 0 | 0 | 0 | 4 | 8 | 10 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.01 | 1118 | 0.190 | 0.55 | 4.264e-07 | 0 | 0 | 0 | 1 | 7 | 9 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.02 | 577 | 0.180 | 0.52 | 2.996e-07 | 0 | 0 | 0 | 0 | 6 | 8 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.04 | 329 | 0.160 | 0.46 | 1.753e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.06 | 258 | 0.140 | 0.40 | 1.690e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.08 | 201 | 0.120 | 0.35 | 1.796e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.1 | 116 | 0.100 | 0.29 | 2.309e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.13 | 137 | 0.070 | 0.20 | 4.432e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.16 | 88 | 0.040 | 0.12 | 5.994e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.18 | 66 | 0.020 | 0.06 | 4.563e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.005 | 3139 | 0.045 | 0.13 | 8.105e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.01 | 1118 | 0.040 | 0.12 | 7.452e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.02 | 577 | 0.030 | 0.09 | 8.201e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.04 | 329 | 0.010 | 0.03 | 2.236e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.005 | 3139 | 0.095 | 0.27 | 1.313e-06 | 0 | 0 | 0 | 1 | 4 | 6 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.01 | 1118 | 0.090 | 0.26 | 9.267e-07 | 0 | 0 | 0 | 0 | 3 | 5 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.02 | 577 | 0.080 | 0.23 | 6.530e-07 | 0 | 0 | 0 | 0 | 1 | 4 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.04 | 329 | 0.060 | 0.17 | 6.318e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.06 | 258 | 0.040 | 0.12 | 1.408e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.08 | 201 | 0.020 | 0.06 | 3.488e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.005 | 3139 | 0.095 | 0.27 | 5.817e-07 | 0 | 0 | 0 | 0 | 4 | 6 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.01 | 1118 | 0.090 | 0.26 | 4.239e-07 | 0 | 0 | 0 | 0 | 1 | 5 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.02 | 577 | 0.080 | 0.23 | 3.275e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.04 | 329 | 0.060 | 0.17 | 3.604e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.06 | 258 | 0.040 | 0.12 | 7.664e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.08 | 201 | 0.020 | 0.06 | 1.682e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.005 | 3139 | 0.195 | 0.56 | 4.777e-07 | 0 | 0 | 0 | 4 | 6 | 8 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.01 | 1118 | 0.190 | 0.55 | 4.085e-07 | 0 | 0 | 0 | 4 | 7 | 8 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.02 | 577 | 0.180 | 0.52 | 3.230e-07 | 0 | 0 | 0 | 3 | 6 | 8 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.04 | 329 | 0.160 | 0.46 | 1.692e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 4 | 7 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.06 | 258 | 0.140 | 0.40 | 1.043e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.08 | 201 | 0.120 | 0.35 | 9.788e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.1 | 116 | 0.100 | 0.29 | 1.077e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.13 | 137 | 0.070 | 0.20 | 2.497e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.16 | 88 | 0.040 | 0.12 | 5.677e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.18 | 66 | 0.020 | 0.06 | 4.282e-06 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.005 | 3139 | 0.195 | 0.56 | 5.313e-07 | 0 | 0 | 0 | 5 | 8 | 9 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.01 | 1118 | 0.190 | 0.55 | 4.273e-07 | 0 | 0 | 0 | 1 | 7 | 9 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.02 | 577 | 0.180 | 0.52 | 2.972e-07 | 0 | 0 | 0 | 0 | 6 | 8 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.04 | 329 | 0.160 | 0.46 | 1.665e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.06 | 258 | 0.140 | 0.40 | 1.560e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.08 | 201 | 0.120 | 0.35 | 1.682e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.1 | 116 | 0.100 | 0.29 | 2.098e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.13 | 137 | 0.070 | 0.20 | 3.975e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.16 | 88 | 0.040 | 0.12 | 5.467e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.18 | 66 | 0.020 | 0.06 | 4.299e-06 | 0 | 0 | 0 | 0 | 0 | 0 |

Differential fit lag = S·(d − r) + b over shells with d − r ≥ 0.08 (prediction S = 2.88 steps per chord, speed 1):

| case | threshold | shells | S | b | R² | implied speed |
|---|---|---|---|---|---|---|
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 2% | 3 | all shells at lag 0 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 5% | 3 | all shells at lag 0 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 10% | 3 | all shells at lag 0 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 25% | 3 | 57 | -4.7 | 0.571 | 0.050 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 50% | 3 | 200 | -15.0 | 1.000 | 0.014 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 75% | 3 | 129 | -6.4 | 0.964 | 0.022 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 2% | 3 | all shells at lag 0 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 5% | 3 | all shells at lag 0 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 10% | 3 | all shells at lag 0 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 25% | 3 | all shells at lag 0 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 50% | 3 | 243 | -19.8 | 0.794 | 0.012 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 75% | 3 | 343 | -26.3 | 0.980 | 0.008 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 2% | 3 | all shells at lag 0 |  |  |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 5% | 3 | all shells at lag 0 |  |  |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 10% | 3 | all shells at lag 0 |  |  |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 25% | 3 | 214 | -37.4 | 0.893 | 0.013 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 50% | 3 | all shells at lag 7 |  |  |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 75% | 3 | 14 | 5.6 | 0.036 | 0.202 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 2% | 5 | all shells at lag 0 |  |  |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 5% | 5 | all shells at lag 0 |  |  |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 10% | 5 | all shells at lag 0 |  |  |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 25% | 5 | 24 | -2.8 | 0.371 | 0.120 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 50% | 5 | 88 | -9.6 | 0.970 | 0.033 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 75% | 5 | 112 | -12.2 | 0.973 | 0.026 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 2% | 3 | all shells at lag 0 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 5% | 3 | all shells at lag 0 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 10% | 3 | all shells at lag 0 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 25% | 3 | 57 | -4.7 | 0.571 | 0.050 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 50% | 3 | 200 | -15.0 | 1.000 | 0.014 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 75% | 3 | 129 | -6.4 | 0.964 | 0.022 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 2% | 3 | all shells at lag 0 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 5% | 3 | all shells at lag 0 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 10% | 3 | all shells at lag 0 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 25% | 3 | all shells at lag 0 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 50% | 3 | 243 | -19.8 | 0.794 | 0.012 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 75% | 3 | 343 | -26.3 | 0.980 | 0.008 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 2% | 3 | all shells at lag 0 |  |  |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 5% | 3 | all shells at lag 0 |  |  |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 10% | 3 | all shells at lag 0 |  |  |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 25% | 3 | 71 | -9.8 | 0.893 | 0.040 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 50% | 3 | 14 | 3.6 | 0.036 | 0.202 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 75% | 3 | all shells at lag 8 |  |  |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 2% | 4 | all shells at lag 0 |  |  |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 5% | 4 | all shells at lag 0 |  |  |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 10% | 4 | all shells at lag 0 |  |  |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 25% | 4 | 28 | -3.2 | 0.277 | 0.103 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 50% | 4 | 80 | -8.1 | 0.991 | 0.036 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 75% | 4 | 98 | -9.7 | 0.997 | 0.030 |

## III.5 Surface pressure and skin friction

RMS and maximum over the 370 airfoil nodes of the difference to the reference at the same step (`rans_surface_norms.csv`). Cf_x is the streamwise wall-shear coefficient. The minimum-Cf_x node of the reference is listed as a separation check (a negative Cf_x at the leading edge is the stagnation region below the nose at AoA 6°, not separation).

| case | when | step | RMS dCf_x | max dCf_x | x of max | RMS dCp | max dCp | x of max | max dy+ | min Cf_x ref | x | max y+ ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | step 9 | 9 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 ConsGalerkinProj Drift | step 200 | 200 | 7.379e-06 | 4.185e-05 | 0.1535 | 1.863e-03 | 8.473e-03 | 0.2284 | 1.820e-02 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 ConsGalerkinProj Drift | last | 399 | 5.621e-06 | 2.160e-05 | 0.9086 | 1.492e-03 | 8.617e-03 | 0.6276 | 8.245e-03 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | step 9 | 9 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | step 200 | 200 | 5.648e-06 | 2.633e-05 | 0.6640 | 1.684e-03 | 9.838e-03 | 0.7316 | 5.543e-03 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | last | 399 | 7.648e-06 | 3.876e-05 | 0.1646 | 1.833e-03 | 9.374e-03 | 0.3233 | 1.493e-02 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 NN Drift | step 9 | 9 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 NN Drift | step 200 | 200 | 7.462e-06 | 4.216e-05 | 0.1535 | 1.839e-03 | 8.130e-03 | 0.2284 | 1.742e-02 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 NN Drift | last | 399 | 5.748e-06 | 2.149e-05 | 0.8919 | 1.544e-03 | 8.971e-03 | 0.6276 | 8.379e-03 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 NN Drift Seed2 | step 9 | 9 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 NN Drift Seed2 | step 200 | 200 | 6.058e-06 | 2.542e-05 | 0.6640 | 1.910e-03 | 9.743e-03 | 0.2284 | 5.869e-03 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.01-20 NN Drift Seed2 | last | 399 | 7.638e-06 | 4.078e-05 | 0.3062 | 1.873e-03 | 1.065e-02 | 0.3233 | 1.421e-02 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | event+1 | 201 | 6.708e-08 | 2.668e-07 | 0.1646 | 1.364e-04 | 6.193e-04 | 0.3957 | 2.850e-04 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | event+50 | 250 | 6.644e-07 | 2.094e-06 | 0.0099 | 1.510e-04 | 5.640e-04 | 0.4533 | 9.990e-04 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | last | 399 | 7.909e-07 | 2.483e-06 | 0.0099 | 1.573e-04 | 5.744e-04 | 0.4533 | 1.189e-03 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | event-1 | 199 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | event+1 | 201 | 2.670e-08 | 1.737e-07 | 0.4533 | 9.072e-05 | 5.053e-04 | 0.3957 | 1.740e-04 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | event+50 | 250 | 6.221e-07 | 3.854e-06 | 0.0087 | 3.267e-04 | 8.401e-04 | 0.0000 | 2.650e-03 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | last | 399 | 7.600e-07 | 5.105e-06 | 0.0087 | 3.656e-04 | 9.356e-04 | 0.0000 | 3.507e-03 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | event+1 | 201 | 1.163e-08 | 7.544e-08 | 0.4533 | 4.832e-05 | 2.783e-04 | 0.3957 | 5.737e-05 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | event+50 | 250 | 2.110e-07 | 6.455e-07 | 0.0087 | 1.301e-04 | 3.204e-04 | 0.0014 | 4.451e-04 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | last | 399 | 3.005e-07 | 1.049e-06 | 0.0087 | 1.567e-04 | 4.005e-04 | 0.0014 | 7.233e-04 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | event-1 | 399 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | event+1 | 401 | 7.377e-10 | 2.561e-09 | 0.9958 | 1.644e-05 | 3.569e-05 | 0.9958 | 4.232e-06 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | event+50 | 450 | 3.761e-07 | 1.481e-06 | 0.0048 | 1.089e-04 | 2.969e-04 | 0.0041 | 3.083e-04 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | last | 799 | 5.700e-07 | 1.585e-06 | 0.0048 | 1.177e-04 | 3.304e-04 | 0.0035 | 6.187e-04 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | event+1 | 201 | 1.245e-09 | 4.424e-09 | 0.9980 | 2.406e-05 | 6.889e-05 | 0.9853 | 2.176e-06 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | event+50 | 250 | 3.817e-07 | 1.517e-06 | 0.0048 | 1.201e-04 | 3.453e-04 | 0.0030 | 5.340e-04 | -2.7840e-03 | 0.0014 | 2.765 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | last | 399 | 4.688e-07 | 1.579e-06 | 0.0048 | 1.315e-04 | 3.695e-04 | 0.0041 | 8.940e-04 | -2.7840e-03 | 0.0014 | 2.765 |
| SST 0.01-20 ConsGalerkinProj Drift | step 9 | 9 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 ConsGalerkinProj Drift | step 200 | 200 | 8.339e-06 | 5.032e-05 | 0.1535 | 2.214e-03 | 1.246e-02 | 0.7769 | 1.570e-02 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 ConsGalerkinProj Drift | last | 399 | 7.046e-06 | 2.379e-05 | 0.4925 | 1.519e-03 | 7.863e-03 | 0.5121 | 1.140e-02 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | step 9 | 9 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | step 200 | 200 | 7.142e-06 | 2.371e-05 | 0.8919 | 1.753e-03 | 8.872e-03 | 0.2284 | 6.445e-03 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | last | 399 | 9.689e-06 | 5.716e-05 | 0.3062 | 2.074e-03 | 1.214e-02 | 0.3233 | 1.427e-02 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 NN Drift | step 9 | 9 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 NN Drift | step 200 | 200 | 8.231e-06 | 4.768e-05 | 0.1535 | 2.346e-03 | 1.200e-02 | 0.7769 | 1.418e-02 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 NN Drift | last | 399 | 7.002e-06 | 2.346e-05 | 0.5121 | 1.510e-03 | 7.654e-03 | 0.5121 | 1.126e-02 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 NN Drift Seed2 | step 9 | 9 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 NN Drift Seed2 | step 200 | 200 | 7.375e-06 | 3.104e-05 | 0.7910 | 1.700e-03 | 8.992e-03 | 0.8173 | 7.366e-03 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.01-20 NN Drift Seed2 | last | 399 | 9.749e-06 | 5.736e-05 | 0.3062 | 2.166e-03 | 1.269e-02 | 0.3233 | 1.462e-02 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | event+1 | 201 | 6.639e-08 | 2.529e-07 | 0.1646 | 1.356e-04 | 6.201e-04 | 0.3957 | 3.375e-04 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | event+50 | 250 | 6.304e-07 | 2.127e-06 | 0.0099 | 1.521e-04 | 5.705e-04 | 0.4533 | 1.012e-03 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | last | 399 | 7.098e-07 | 2.486e-06 | 0.0099 | 1.639e-04 | 5.828e-04 | 0.4533 | 1.186e-03 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | event-1 | 199 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | event+1 | 201 | 2.818e-08 | 1.814e-07 | 0.4533 | 8.993e-05 | 5.077e-04 | 0.3957 | 1.622e-04 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | event+50 | 250 | 6.164e-07 | 4.064e-06 | 0.0087 | 3.248e-04 | 8.367e-04 | 0.0000 | 2.822e-03 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | last | 399 | 7.733e-07 | 5.181e-06 | 0.0087 | 3.612e-04 | 9.255e-04 | 0.0000 | 3.594e-03 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | event+1 | 201 | 1.296e-08 | 8.196e-08 | 0.4533 | 4.903e-05 | 2.769e-04 | 0.3957 | 6.092e-05 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | event+50 | 250 | 1.880e-07 | 6.866e-07 | 0.0087 | 1.311e-04 | 3.228e-04 | 0.0014 | 4.781e-04 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | last | 399 | 2.836e-07 | 1.089e-06 | 0.0087 | 1.586e-04 | 4.051e-04 | 0.0014 | 7.581e-04 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | event-1 | 399 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | event+1 | 401 | 1.004e-09 | 6.519e-09 | 0.0009 | 1.557e-05 | 3.089e-05 | 0.9909 | 3.606e-06 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | event+50 | 450 | 3.590e-07 | 1.376e-06 | 0.0048 | 1.086e-04 | 3.009e-04 | 0.0041 | 3.416e-04 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | last | 799 | 5.218e-07 | 1.445e-06 | 0.0048 | 1.151e-04 | 3.171e-04 | 0.0041 | 6.434e-04 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | 0.000e+00 | 1.0000 | 0.000e+00 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | event+1 | 201 | 1.321e-09 | 3.842e-09 | 0.9934 | 2.144e-05 | 5.662e-05 | 0.9980 | 4.023e-06 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | event+50 | 250 | 3.769e-07 | 1.454e-06 | 0.0041 | 1.201e-04 | 3.468e-04 | 0.0030 | 5.975e-04 | -2.7153e-03 | 0.0014 | 2.698 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | last | 399 | 4.447e-07 | 1.480e-06 | 0.0048 | 1.304e-04 | 3.648e-04 | 0.0041 | 9.567e-04 | -2.7153e-03 | 0.0014 | 2.698 |

Distributions (Cf_x and Cp of case and reference at the key steps) are in `rans_surface_profiles.csv` and plotted in `Figures/RANS/<case>_Cf_Cp.png`.

## III.6 Turbulence state

**Transfer defects of the transported variables.** For every transfer, the domain integral of the rho-weighted variable on the source mesh and on the target mesh (`rans_transfer_integrals.csv`, from the interpolation's own `Integrals_<step>.csv`). Only the newest state of each event is summarised (the BDF2 level is transferred too and behaves the same). The last column is the change of the source integral between consecutive events: it is what the transport equations did in the 10 steps in between, and it is listed so that it is not mistaken for a transfer defect.

| case | field | events | max |defect| | median |defect| | min change between events | max change between events | net change first→last event |
|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | Density | 40 | 4.4e-13 | 2.5e-13 | -3.596e-11 | 3.783e-11 | 5.229e-11 |
| SA 0.01-20 ConsGalerkinProj Drift | Energy | 40 | 5.8e-13 | 3.0e-13 | -4.942e-11 | 5.521e-11 | 9.900e-11 |
| SA 0.01-20 ConsGalerkinProj Drift | Nu_Tilde | 40 | 6.6e-09 | 2.3e-09 | -6.121e-07 | 1.305e-07 | -7.299e-06 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | Density | 40 | 4.5e-13 | 2.3e-13 | -2.337e-11 | 3.339e-11 | 3.380e-11 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | Energy | 40 | 5.4e-13 | 3.0e-13 | -3.058e-11 | 4.771e-11 | 7.763e-11 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | Nu_Tilde | 40 | 6.9e-09 | 1.9e-09 | -4.511e-07 | 7.276e-07 | -3.376e-06 |
| SA 0.01-20 NN Drift | Density | 40 | 6.8e-11 | 2.0e-11 | -7.130e-11 | 8.913e-11 | 2.257e-10 |
| SA 0.01-20 NN Drift | Energy | 40 | 5.7e-11 | 1.7e-11 | -5.416e-11 | 1.189e-10 | 5.573e-10 |
| SA 0.01-20 NN Drift | Nu_Tilde | 40 | 1.7e-06 | 5.3e-07 | -1.695e-06 | 1.271e-06 | -1.786e-06 |
| SA 0.01-20 NN Drift Seed2 | Density | 40 | 8.4e-11 | 2.0e-11 | -8.709e-11 | 8.208e-11 | -1.476e-10 |
| SA 0.01-20 NN Drift Seed2 | Energy | 40 | 9.5e-11 | 2.2e-11 | -9.967e-11 | 8.306e-11 | -1.004e-10 |
| SA 0.01-20 NN Drift Seed2 | Nu_Tilde | 40 | 1.5e-06 | 3.5e-07 | -1.649e-06 | 1.621e-06 | -4.631e-06 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | Density | 2 | 2.3e-13 | 1.2e-13 | 3.671e-13 | 3.671e-13 | 3.189e-13 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | Energy | 2 | 3.0e-13 | 1.6e-13 | 2.428e-13 | 2.428e-13 | 1.781e-13 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | Nu_Tilde | 2 | 7.5e-12 | 5.5e-12 | 3.735e-07 | 3.735e-07 | 3.735e-07 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | Density | 2 | 9.9e-14 | 6.2e-14 | 2.298e-13 | 2.298e-13 | 3.189e-13 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | Energy | 2 | 1.3e-13 | 7.2e-14 | -2.675e-13 | -2.675e-13 | -3.560e-13 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | Nu_Tilde | 2 | 1.3e-12 | 7.4e-13 | 4.801e-07 | 4.801e-07 | 4.801e-07 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | Density | 2 | 1.2e-14 | 1.1e-14 | 3.036e-13 | 3.036e-13 | 3.189e-13 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | Energy | 2 | 1.5e-14 | 7.8e-15 | -5.023e-13 | -5.023e-13 | -5.341e-13 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | Nu_Tilde | 2 | 2.5e-13 | 1.6e-13 | -3.638e-08 | -3.638e-08 | -3.638e-08 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | Density | 2 | 9.3e-16 | 7.5e-16 | -5.174e-12 | -5.174e-12 | -5.101e-12 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | Energy | 2 | 1.2e-15 | 7.8e-16 | -5.794e-12 | -5.794e-12 | -5.876e-12 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | Nu_Tilde | 2 | 5.9e-15 | 3.0e-15 | -2.599e-07 | -2.599e-07 | -2.599e-07 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | Density | 2 | 4.4e-15 | 2.4e-15 | -7.258e-13 | -7.258e-13 | -6.377e-13 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | Energy | 2 | 1.2e-15 | 6.5e-16 | 1.474e-12 | 1.474e-12 | 1.424e-12 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | Nu_Tilde | 2 | 4.6e-15 | 4.4e-15 | -8.473e-08 | -8.473e-08 | -8.473e-08 |
| SST 0.01-20 ConsGalerkinProj Drift | Density | 40 | 4.5e-13 | 2.6e-13 | -3.609e-11 | 3.601e-11 | 6.122e-11 |
| SST 0.01-20 ConsGalerkinProj Drift | Energy | 40 | 5.7e-13 | 3.1e-13 | -4.976e-11 | 5.214e-11 | 9.704e-11 |
| SST 0.01-20 ConsGalerkinProj Drift | Turb_Kin_Energy | 40 | 4.2e-08 | 2.0e-08 | -3.764e-06 | 1.506e-06 | -2.603e-05 |
| SST 0.01-20 ConsGalerkinProj Drift | Omega | 40 | 2.1e-09 | 5.1e-10 | -1.080e-07 | 1.833e-07 | 5.676e-07 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Density | 40 | 4.5e-13 | 2.4e-13 | -2.477e-11 | 3.260e-11 | 5.452e-11 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Energy | 40 | 5.6e-13 | 3.1e-13 | -3.295e-11 | 4.614e-11 | 9.330e-11 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Turb_Kin_Energy | 40 | 4.8e-08 | 1.8e-08 | -7.639e-06 | 2.077e-05 | -1.496e-05 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Omega | 40 | 1.6e-09 | 5.7e-10 | -5.905e-07 | 1.161e-06 | 8.980e-07 |
| SST 0.01-20 NN Drift | Density | 40 | 6.6e-11 | 2.0e-11 | -6.935e-11 | 9.011e-11 | 2.238e-10 |
| SST 0.01-20 NN Drift | Energy | 40 | 5.7e-11 | 1.7e-11 | -5.314e-11 | 1.154e-10 | 5.329e-10 |
| SST 0.01-20 NN Drift | Turb_Kin_Energy | 40 | 2.2e-06 | 5.5e-07 | -4.008e-06 | 2.643e-06 | -2.401e-05 |
| SST 0.01-20 NN Drift | Omega | 40 | 7.9e-08 | 2.1e-08 | -1.405e-07 | 1.925e-07 | 6.350e-07 |
| SST 0.01-20 NN Drift Seed2 | Density | 40 | 8.3e-11 | 1.9e-11 | -8.576e-11 | 7.823e-11 | -1.253e-10 |
| SST 0.01-20 NN Drift Seed2 | Energy | 40 | 9.4e-11 | 2.2e-11 | -9.870e-11 | 8.767e-11 | -8.654e-11 |
| SST 0.01-20 NN Drift Seed2 | Turb_Kin_Energy | 40 | 1.8e-06 | 4.1e-07 | -6.824e-06 | 1.915e-05 | -2.447e-05 |
| SST 0.01-20 NN Drift Seed2 | Omega | 40 | 6.0e-08 | 2.0e-08 | -5.808e-07 | 1.049e-06 | 6.371e-07 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Density | 2 | 2.3e-13 | 1.2e-13 | 2.166e-12 | 2.166e-12 | 2.232e-12 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Energy | 2 | 2.9e-13 | 1.5e-13 | 9.687e-13 | 9.687e-13 | 8.902e-13 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Turb_Kin_Energy | 2 | 6.6e-12 | 5.2e-12 | -1.034e-06 | -1.034e-06 | -1.034e-06 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Omega | 2 | 1.9e-11 | 1.2e-11 | -3.163e-08 | -3.163e-08 | -3.163e-08 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Density | 2 | 9.8e-14 | 6.0e-14 | 1.574e-12 | 1.574e-12 | 1.594e-12 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Energy | 2 | 1.3e-13 | 7.0e-14 | -8.473e-14 | -8.473e-14 | -1.781e-13 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Turb_Kin_Energy | 2 | 2.2e-11 | 1.1e-11 | 2.648e-06 | 2.648e-06 | 2.648e-06 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Omega | 2 | 1.6e-12 | 1.5e-12 | 7.905e-08 | 7.905e-08 | 7.906e-08 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Density | 2 | 1.2e-14 | 1.1e-14 | 1.456e-12 | 1.456e-12 | 1.275e-12 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Energy | 2 | 1.5e-14 | 8.0e-15 | -3.563e-13 | -3.563e-13 | -3.562e-13 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Turb_Kin_Energy | 2 | 7.6e-14 | 6.4e-14 | 5.143e-07 | 5.143e-07 | 5.143e-07 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Omega | 2 | 6.2e-13 | 3.6e-13 | 2.238e-08 | 2.238e-08 | 2.238e-08 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Density | 2 | 2.1e-15 | 1.7e-15 | -1.636e-12 | -1.636e-12 | -1.594e-12 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Energy | 2 | 2.6e-16 | 1.3e-16 | -4.904e-12 | -4.904e-12 | -4.986e-12 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Turb_Kin_Energy | 2 | 4.5e-14 | 4.5e-14 | -2.820e-06 | -2.820e-06 | -2.820e-06 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Omega | 2 | 1.4e-13 | 1.4e-13 | -5.957e-08 | -5.957e-08 | -5.957e-08 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Density | 2 | 2.8e-15 | 1.4e-15 | 5.946e-13 | 5.946e-13 | 6.377e-13 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Energy | 2 | 1.4e-15 | 8.4e-16 | 1.618e-12 | 1.618e-12 | 1.602e-12 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Turb_Kin_Energy | 2 | 4.4e-10 | 2.2e-10 | 7.013e-07 | 7.013e-07 | 7.013e-07 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Omega | 2 | 1.4e-13 | 1.4e-13 | 1.316e-08 | 1.316e-08 | 1.316e-08 |

**Minimum of the transferred turbulence variables** over all transfers of each case (source mesh and target mesh, from the interpolation's `Minimums_<step>.csv`). These are the rho-weighted values the interpolation checks (rho*nu_tilde, rho*k, rho*omega), so the free-stream floor of k appears as 0.98e-10 rather than 1e-10; the clipping of III.6 acts on the primitive variables after the division by rho, before this check. "files with target < 0" counts the transfers whose interpolated field contained a negative value before that clipping (SST) or, for SA, which runs without clipping, in the returned state:

| case | field | transfers | min source | min target | files with target < 0 |
|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | Nu_Tilde | 80 | 0.000e+00 | -6.150e-17 | 70 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | Nu_Tilde | 80 | 0.000e+00 | -6.377e-17 | 62 |
| SA 0.01-20 NN Drift | Nu_Tilde | 80 | 0.000e+00 | -2.954e-08 | 5 |
| SA 0.01-20 NN Drift Seed2 | Nu_Tilde | 80 | 0.000e+00 | -2.698e-09 | 2 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | Nu_Tilde | 4 | 0.000e+00 | -4.980e-18 | 4 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | Nu_Tilde | 4 | 0.000e+00 | -5.372e-18 | 4 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | Nu_Tilde | 4 | 0.000e+00 | -4.915e-18 | 4 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | Nu_Tilde | 4 | 0.000e+00 | -4.922e-18 | 4 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | Nu_Tilde | 4 | 0.000e+00 | -4.792e-18 | 4 |
| SST 0.01-20 ConsGalerkinProj Drift | Turb_Kin_Energy | 80 | 9.814e-11 | 9.814e-11 | 0 |
| SST 0.01-20 ConsGalerkinProj Drift | Omega | 80 | 7.660e-03 | 3.380e-03 | 0 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Turb_Kin_Energy | 80 | 9.814e-11 | 9.814e-11 | 0 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Omega | 80 | 7.660e-03 | 9.994e-05 | 0 |
| SST 0.01-20 NN Drift | Turb_Kin_Energy | 80 | 9.814e-11 | 9.814e-11 | 0 |
| SST 0.01-20 NN Drift | Omega | 80 | 7.660e-03 | 7.660e-03 | 0 |
| SST 0.01-20 NN Drift Seed2 | Turb_Kin_Energy | 80 | 9.814e-11 | 9.814e-11 | 0 |
| SST 0.01-20 NN Drift Seed2 | Omega | 80 | 7.660e-03 | 1.000e-04 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Turb_Kin_Energy | 4 | 9.814e-11 | 9.814e-11 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Omega | 4 | 7.660e-03 | 7.660e-03 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Turb_Kin_Energy | 4 | 9.814e-11 | 9.814e-11 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Omega | 4 | 7.660e-03 | 7.660e-03 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Turb_Kin_Energy | 4 | 9.814e-11 | 9.814e-11 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Omega | 4 | 7.660e-03 | 7.660e-03 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Turb_Kin_Energy | 4 | 9.814e-11 | 9.814e-11 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Omega | 4 | 7.660e-03 | 7.660e-03 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Turb_Kin_Energy | 4 | 9.814e-11 | 9.814e-11 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Omega | 4 | 7.660e-03 | 7.660e-03 | 0 |

**Clipping to the solver's lower limits inside the interpolation** (`clipFieldsToMinimum` / `FieldsMinimum` in InterpSolution_Params.py): after each transfer, values of the listed fields below their floor are raised to it. The floors are the ones SU2 applies to its SST variables (k ≥ 1e-10, ω ≥ 1e-4). Per case and field: number of transfers, transfers with at least one raised node, largest number of raised nodes in one transfer, lowest value seen before clipping, largest total amount added in one transfer. The SA cases ran with the option off (their lowest transferred nu_tilde was of order -1e-16, see the minima table above).

| case | setting | field | transfers | floor | with raises | max nodes | lowest before | max added |
|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | [] | - | 80 | - | - | - | - | - |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | [] | - | 80 | - | - | - | - | - |
| SA 0.01-20 NN Drift | [] | - | 80 | - | - | - | - | - |
| SA 0.01-20 NN Drift Seed2 | [] | - | 80 | - | - | - | - | - |
| SA 0.05-0.10 ConsGalerkinProj OneShot | [] | - | 4 | - | - | - | - | - |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | [] | - | 4 | - | - | - | - | - |
| SA 0.10-0.15 ConsGalerkinProj OneShot | [] | - | 4 | - | - | - | - | - |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | [] | - | 4 | - | - | - | - | - |
| SA 0.20-0.25 ConsGalerkinProj OneShot | [] | - | 4 | - | - | - | - | - |
| SST 0.01-20 ConsGalerkinProj Drift | ['Turb_Kin_Energy', 'Omega'] | Omega | 80 | 0.0001 | 0 | 0 | 3.38e-03 | 0.00e+00 |
| SST 0.01-20 ConsGalerkinProj Drift | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 80 | 1e-10 | 80 | 191 | -2.21e-07 | 5.24e-07 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | ['Turb_Kin_Energy', 'Omega'] | Omega | 80 | 0.0001 | 2 | 2 | -7.83e-02 | 9.61e-02 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 80 | 1e-10 | 80 | 337 | -4.02e-07 | 1.32e-06 |
| SST 0.01-20 NN Drift | ['Turb_Kin_Energy', 'Omega'] | Omega | 80 | 0.0001 | 0 | 0 | 7.66e-03 | 0.00e+00 |
| SST 0.01-20 NN Drift | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 80 | 1e-10 | 80 | 13 | -1.98e-07 | 3.71e-07 |
| SST 0.01-20 NN Drift Seed2 | ['Turb_Kin_Energy', 'Omega'] | Omega | 80 | 0.0001 | 2 | 1 | -3.22e-02 | 3.23e-02 |
| SST 0.01-20 NN Drift Seed2 | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 80 | 1e-10 | 80 | 20 | -3.49e-07 | 9.69e-07 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Omega | 4 | 0.0001 | 0 | 0 | 7.66e-03 | 0.00e+00 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 4 | 1e-10 | 4 | 168 | 1.00e-10 | 4.51e-16 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | ['Turb_Kin_Energy', 'Omega'] | Omega | 4 | 0.0001 | 0 | 0 | 7.66e-03 | 0.00e+00 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 4 | 1e-10 | 4 | 169 | -2.06e-10 | 3.06e-10 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Omega | 4 | 0.0001 | 0 | 0 | 7.66e-03 | 0.00e+00 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 4 | 1e-10 | 4 | 168 | 1.00e-10 | 4.50e-16 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | ['Turb_Kin_Energy', 'Omega'] | Omega | 4 | 0.0001 | 0 | 0 | 7.66e-03 | 0.00e+00 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 4 | 1e-10 | 4 | 168 | 1.00e-10 | 4.50e-16 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Omega | 4 | 0.0001 | 0 | 0 | 7.66e-03 | 0.00e+00 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 4 | 1e-10 | 4 | 169 | -2.21e-09 | 2.31e-09 |

**SA: node-wise discrepancy on never-displaced nodes** (`rans_turb_state_norms_SA.csv`). Nodes whose coordinates coincide with the reference to 1e-12 are compared directly; "band" = wall distance < 0.01 (the wall-layer nodes), "outer" = the other fixed nodes. Values are RMS and max of |case − reference| divided by the reference maximum over the same node set.

| case | step | fixed nodes | of which band | moved nodes | Eddy_Viscosity band RMS | band max | outer RMS | outer max | Nu_Tilde band RMS | band max | outer RMS | outer max |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | 10 | 10669 | 8590 | 12602 | 1.345e-04 | 3.165e-03 | 6.930e-05 | 1.727e-03 | 1.334e-04 | 3.182e-03 | 2.857e-05 | 7.290e-04 |
| SA 0.01-20 ConsGalerkinProj Drift | 200 | 10669 | 8590 | 12602 | 7.702e-04 | 1.074e-02 | 2.123e-04 | 6.027e-03 | 7.694e-04 | 1.073e-02 | 9.747e-05 | 2.832e-03 |
| SA 0.01-20 ConsGalerkinProj Drift | 395 | 10669 | 8590 | 12602 | 1.165e-03 | 1.055e-02 | 3.043e-04 | 6.372e-03 | 1.165e-03 | 1.055e-02 | 1.352e-04 | 2.893e-03 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | 10 | 10669 | 8590 | 12602 | 1.680e-04 | 3.917e-03 | 8.684e-05 | 2.548e-03 | 1.673e-04 | 3.918e-03 | 3.398e-05 | 7.917e-04 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | 200 | 10669 | 8590 | 12602 | 1.066e-03 | 1.229e-02 | 3.085e-04 | 8.039e-03 | 1.066e-03 | 1.229e-02 | 1.243e-04 | 2.456e-03 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | 395 | 10669 | 8590 | 12602 | 1.416e-03 | 1.635e-02 | 2.830e-04 | 5.755e-03 | 1.416e-03 | 1.636e-02 | 1.101e-04 | 1.755e-03 |
| SA 0.01-20 NN Drift | 10 | 10669 | 8590 | 12602 | 1.265e-04 | 2.855e-03 | 4.289e-06 | 1.479e-04 | 1.256e-04 | 2.863e-03 | 1.903e-06 | 4.988e-05 |
| SA 0.01-20 NN Drift | 200 | 10669 | 8590 | 12602 | 6.811e-04 | 9.931e-03 | 3.004e-04 | 6.959e-03 | 6.804e-04 | 9.928e-03 | 1.248e-04 | 2.749e-03 |
| SA 0.01-20 NN Drift | 395 | 10669 | 8590 | 12602 | 1.036e-03 | 1.045e-02 | 3.536e-04 | 8.128e-03 | 1.036e-03 | 1.045e-02 | 1.426e-04 | 2.469e-03 |
| SA 0.01-20 NN Drift Seed2 | 10 | 10669 | 8590 | 12602 | 1.672e-04 | 3.754e-03 | 6.745e-06 | 2.726e-04 | 1.666e-04 | 3.756e-03 | 2.471e-06 | 8.295e-05 |
| SA 0.01-20 NN Drift Seed2 | 200 | 10669 | 8590 | 12602 | 1.144e-03 | 1.351e-02 | 2.521e-04 | 5.863e-03 | 1.145e-03 | 1.351e-02 | 1.077e-04 | 2.754e-03 |
| SA 0.01-20 NN Drift Seed2 | 395 | 10669 | 8590 | 12602 | 8.336e-04 | 1.048e-02 | 2.263e-04 | 4.852e-03 | 8.332e-04 | 1.048e-02 | 9.684e-05 | 2.203e-03 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 200 | 22047 | 8590 | 1224 | 1.825e-06 | 2.785e-05 | 1.141e-04 | 7.982e-03 | 1.767e-06 | 2.808e-05 | 1.139e-04 | 7.995e-03 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 250 | 22047 | 8590 | 1224 | 5.488e-05 | 3.274e-04 | 1.945e-04 | 1.094e-02 | 5.482e-05 | 3.262e-04 | 1.935e-04 | 1.096e-02 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 300 | 22047 | 8590 | 1224 | 5.850e-05 | 3.549e-04 | 1.949e-04 | 1.098e-02 | 5.835e-05 | 3.542e-04 | 1.938e-04 | 1.100e-02 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 395 | 22047 | 8590 | 1224 | 5.365e-05 | 3.230e-04 | 1.942e-04 | 1.094e-02 | 5.346e-05 | 3.224e-04 | 1.932e-04 | 1.096e-02 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 200 | 22508 | 8590 | 763 | 9.909e-07 | 1.561e-05 | 1.193e-04 | 1.130e-02 | 9.633e-07 | 1.549e-05 | 1.186e-04 | 1.131e-02 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 250 | 22508 | 8590 | 763 | 4.007e-05 | 3.098e-04 | 2.689e-04 | 1.700e-02 | 3.999e-05 | 3.096e-04 | 2.644e-04 | 1.700e-02 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 300 | 22508 | 8590 | 763 | 6.125e-05 | 4.064e-04 | 2.694e-04 | 1.702e-02 | 6.117e-05 | 4.060e-04 | 2.649e-04 | 1.702e-02 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 395 | 22508 | 8590 | 763 | 7.054e-05 | 4.583e-04 | 2.696e-04 | 1.703e-02 | 7.043e-05 | 4.576e-04 | 2.650e-04 | 1.704e-02 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 200 | 22508 | 8590 | 763 | 3.525e-07 | 4.016e-06 | 5.737e-05 | 5.577e-03 | 3.415e-07 | 4.072e-06 | 5.706e-05 | 5.584e-03 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 250 | 22508 | 8590 | 763 | 1.046e-05 | 6.028e-05 | 1.233e-04 | 8.118e-03 | 1.045e-05 | 6.027e-05 | 1.211e-04 | 8.123e-03 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 300 | 22508 | 8590 | 763 | 1.994e-05 | 1.163e-04 | 1.236e-04 | 8.133e-03 | 1.993e-05 | 1.162e-04 | 1.215e-04 | 8.137e-03 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 395 | 22508 | 8590 | 763 | 3.051e-05 | 1.720e-04 | 1.239e-04 | 8.132e-03 | 3.048e-05 | 1.718e-04 | 1.217e-04 | 8.136e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 400 | 22877 | 8590 | 394 | 3.317e-08 | 5.126e-07 | 3.997e-05 | 3.867e-03 | 3.033e-08 | 4.659e-07 | 3.971e-05 | 3.855e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 450 | 22877 | 8590 | 394 | 3.700e-05 | 2.707e-04 | 9.639e-05 | 6.096e-03 | 3.703e-05 | 2.712e-04 | 9.501e-05 | 6.076e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 500 | 22877 | 8590 | 394 | 5.263e-05 | 3.317e-04 | 9.847e-05 | 6.115e-03 | 5.263e-05 | 3.318e-04 | 9.710e-05 | 6.095e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 795 | 22877 | 8590 | 394 | 4.991e-05 | 2.368e-04 | 9.296e-05 | 6.079e-03 | 4.985e-05 | 2.367e-04 | 9.151e-05 | 6.059e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 200 | 22877 | 8590 | 394 | 1.408e-07 | 2.022e-06 | 4.644e-05 | 3.478e-03 | 1.562e-07 | 2.107e-06 | 4.612e-05 | 3.476e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 250 | 22877 | 8590 | 394 | 3.052e-05 | 2.542e-04 | 9.714e-05 | 5.074e-03 | 3.060e-05 | 2.548e-04 | 9.621e-05 | 5.071e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 300 | 22877 | 8590 | 394 | 4.957e-05 | 3.387e-04 | 1.020e-04 | 5.126e-03 | 4.963e-05 | 3.392e-04 | 1.012e-04 | 5.123e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 395 | 22877 | 8590 | 394 | 5.514e-05 | 3.484e-04 | 1.006e-04 | 5.150e-03 | 5.513e-05 | 3.484e-04 | 9.968e-05 | 5.147e-03 |

**SST: node-wise discrepancy on never-displaced nodes** (`rans_turb_state_norms_SST.csv`). Nodes whose coordinates coincide with the reference to 1e-12 are compared directly; "band" = wall distance < 0.01 (the wall-layer nodes), "outer" = the other fixed nodes. Values are RMS and max of |case − reference| divided by the reference maximum over the same node set.

| case | step | fixed nodes | of which band | moved nodes | Eddy_Viscosity band RMS | band max | outer RMS | outer max | Turb_Kin_Energy band RMS | band max | outer RMS | outer max | Omega band RMS | band max | outer RMS | outer max |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SST 0.01-20 ConsGalerkinProj Drift | 10 | 10669 | 8590 | 12602 | 1.556e-03 | 4.732e-02 | 6.946e-05 | 1.923e-03 | 1.305e-04 | 2.385e-03 | 1.841e-06 | 4.715e-05 | 2.131e-06 | 1.534e-04 | 7.007e-08 | 1.479e-06 |
| SST 0.01-20 ConsGalerkinProj Drift | 200 | 10669 | 8590 | 12602 | 3.444e-03 | 1.152e-01 | 2.425e-04 | 6.707e-03 | 1.005e-03 | 1.019e-02 | 6.320e-06 | 1.857e-04 | 3.944e-06 | 7.403e-05 | 2.370e-07 | 5.973e-06 |
| SST 0.01-20 ConsGalerkinProj Drift | 395 | 10669 | 8590 | 12602 | 3.042e-03 | 6.666e-02 | 3.378e-04 | 8.758e-03 | 5.951e-04 | 5.818e-03 | 8.572e-06 | 1.896e-04 | 3.226e-06 | 5.969e-05 | 2.907e-07 | 5.615e-06 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | 10 | 10669 | 8590 | 12602 | 1.991e-03 | 5.633e-02 | 8.150e-05 | 2.161e-03 | 1.520e-04 | 3.106e-03 | 2.182e-06 | 5.304e-05 | 1.097e-06 | 3.388e-05 | 8.383e-08 | 1.687e-06 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | 200 | 10669 | 8590 | 12602 | 3.064e-03 | 7.830e-02 | 2.895e-04 | 5.748e-03 | 1.001e-03 | 1.053e-02 | 7.936e-06 | 1.571e-04 | 3.946e-06 | 8.335e-05 | 2.762e-07 | 4.895e-06 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | 395 | 10669 | 8590 | 12602 | 3.601e-03 | 9.260e-02 | 2.600e-04 | 5.262e-03 | 1.472e-03 | 1.728e-02 | 6.793e-06 | 1.095e-04 | 4.080e-06 | 8.606e-05 | 2.343e-07 | 3.380e-06 |
| SST 0.01-20 NN Drift | 10 | 10669 | 8590 | 12602 | 1.535e-03 | 4.776e-02 | 5.740e-06 | 1.939e-04 | 1.262e-04 | 2.082e-03 | 1.285e-07 | 3.690e-06 | 2.081e-06 | 1.496e-04 | 6.819e-09 | 1.123e-07 |
| SST 0.01-20 NN Drift | 200 | 10669 | 8590 | 12602 | 3.447e-03 | 1.206e-01 | 2.984e-04 | 7.678e-03 | 9.197e-04 | 9.327e-03 | 8.146e-06 | 1.805e-04 | 3.910e-06 | 7.472e-05 | 2.821e-07 | 5.717e-06 |
| SST 0.01-20 NN Drift | 395 | 10669 | 8590 | 12602 | 3.069e-03 | 6.678e-02 | 3.233e-04 | 7.451e-03 | 6.055e-04 | 5.533e-03 | 8.567e-06 | 1.595e-04 | 3.251e-06 | 6.076e-05 | 2.930e-07 | 4.765e-06 |
| SST 0.01-20 NN Drift Seed2 | 10 | 10669 | 8590 | 12602 | 1.995e-03 | 5.645e-02 | 5.479e-06 | 1.342e-04 | 1.497e-04 | 3.162e-03 | 1.643e-07 | 5.692e-06 | 1.086e-06 | 3.299e-05 | 7.367e-09 | 1.983e-07 |
| SST 0.01-20 NN Drift Seed2 | 200 | 10669 | 8590 | 12602 | 4.385e-03 | 9.574e-02 | 2.623e-04 | 6.444e-03 | 1.256e-03 | 1.062e-02 | 7.039e-06 | 1.796e-04 | 4.100e-06 | 1.162e-04 | 2.496e-07 | 5.892e-06 |
| SST 0.01-20 NN Drift Seed2 | 395 | 10669 | 8590 | 12602 | 2.784e-03 | 8.912e-02 | 2.370e-04 | 6.175e-03 | 1.137e-03 | 9.777e-03 | 6.226e-06 | 1.521e-04 | 4.037e-06 | 8.473e-05 | 2.113e-07 | 4.950e-06 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 200 | 22047 | 8590 | 1224 | 6.308e-06 | 2.053e-04 | 2.035e-04 | 1.680e-02 | 2.750e-06 | 2.729e-05 | 7.836e-05 | 4.692e-03 | 1.210e-07 | 4.634e-06 | 4.878e-05 | 3.511e-03 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 250 | 22047 | 8590 | 1224 | 5.163e-05 | 3.711e-04 | 4.782e-04 | 2.993e-02 | 6.587e-05 | 4.519e-04 | 1.261e-04 | 5.896e-03 | 7.648e-07 | 1.023e-05 | 8.264e-05 | 4.622e-03 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 300 | 22047 | 8590 | 1224 | 6.033e-05 | 3.593e-04 | 4.777e-04 | 2.996e-02 | 6.975e-05 | 4.813e-04 | 1.266e-04 | 5.951e-03 | 7.545e-07 | 1.012e-05 | 8.233e-05 | 4.630e-03 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 395 | 22047 | 8590 | 1224 | 5.976e-05 | 3.538e-04 | 4.774e-04 | 2.997e-02 | 7.208e-05 | 4.938e-04 | 1.269e-04 | 5.950e-03 | 7.212e-07 | 9.427e-06 | 8.253e-05 | 4.632e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 200 | 22508 | 8590 | 763 | 2.867e-06 | 1.032e-04 | 2.548e-04 | 1.842e-02 | 6.858e-07 | 9.756e-06 | 7.680e-05 | 7.446e-03 | 7.833e-08 | 3.348e-06 | 4.577e-05 | 3.784e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 250 | 22508 | 8590 | 763 | 2.368e-05 | 2.850e-04 | 9.463e-04 | 4.927e-02 | 5.892e-05 | 7.725e-04 | 1.555e-04 | 1.104e-02 | 8.113e-07 | 1.023e-05 | 8.794e-05 | 5.601e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 300 | 22508 | 8590 | 763 | 3.302e-05 | 3.261e-04 | 9.476e-04 | 4.926e-02 | 6.472e-05 | 8.606e-04 | 1.555e-04 | 1.103e-02 | 8.256e-07 | 1.023e-05 | 8.785e-05 | 5.602e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 395 | 22508 | 8590 | 763 | 4.937e-05 | 3.487e-04 | 9.494e-04 | 4.928e-02 | 7.241e-05 | 9.456e-04 | 1.555e-04 | 1.104e-02 | 8.153e-07 | 9.664e-06 | 8.761e-05 | 5.604e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 200 | 22508 | 8590 | 763 | 1.197e-06 | 4.463e-05 | 1.278e-04 | 9.229e-03 | 3.669e-07 | 4.630e-06 | 3.711e-05 | 3.666e-03 | 4.805e-08 | 1.466e-06 | 2.225e-05 | 1.868e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 250 | 22508 | 8590 | 763 | 6.396e-06 | 8.661e-05 | 4.686e-04 | 2.561e-02 | 2.792e-05 | 2.564e-04 | 7.167e-05 | 5.207e-03 | 5.112e-07 | 7.242e-06 | 4.007e-05 | 2.669e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 300 | 22508 | 8590 | 763 | 9.416e-06 | 7.294e-05 | 4.692e-04 | 2.561e-02 | 3.097e-05 | 3.125e-04 | 7.192e-05 | 5.207e-03 | 5.178e-07 | 7.242e-06 | 4.037e-05 | 2.668e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 395 | 22508 | 8590 | 763 | 1.697e-05 | 8.998e-05 | 4.698e-04 | 2.561e-02 | 3.469e-05 | 3.641e-04 | 7.197e-05 | 5.206e-03 | 5.046e-07 | 6.783e-06 | 4.023e-05 | 2.667e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 400 | 22877 | 8590 | 394 | 1.935e-07 | 3.858e-06 | 6.748e-05 | 3.286e-03 | 6.557e-08 | 1.449e-06 | 2.039e-05 | 1.984e-03 | 1.361e-08 | 2.299e-07 | 1.233e-05 | 1.178e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 450 | 22877 | 8590 | 394 | 2.320e-05 | 1.505e-04 | 2.297e-04 | 1.427e-02 | 4.699e-05 | 1.949e-04 | 5.881e-05 | 3.053e-03 | 4.671e-07 | 6.323e-06 | 2.892e-05 | 1.836e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 500 | 22877 | 8590 | 394 | 3.785e-05 | 2.465e-04 | 2.332e-04 | 1.429e-02 | 5.157e-05 | 1.886e-04 | 6.078e-05 | 3.053e-03 | 4.708e-07 | 6.438e-06 | 2.908e-05 | 1.834e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 795 | 22877 | 8590 | 394 | 4.939e-05 | 2.425e-04 | 2.326e-04 | 1.431e-02 | 6.002e-05 | 1.736e-04 | 5.191e-05 | 3.042e-03 | 4.257e-07 | 5.518e-06 | 2.730e-05 | 1.830e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 200 | 22877 | 8590 | 394 | 3.793e-07 | 4.549e-06 | 1.097e-04 | 8.726e-03 | 1.180e-07 | 7.953e-07 | 2.324e-05 | 1.763e-03 | 2.649e-08 | 3.449e-07 | 1.443e-05 | 1.135e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 250 | 22877 | 8590 | 394 | 1.878e-05 | 1.468e-04 | 3.494e-04 | 2.712e-02 | 4.648e-05 | 2.041e-04 | 6.055e-05 | 2.555e-03 | 4.882e-07 | 6.668e-06 | 3.117e-05 | 1.682e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 300 | 22877 | 8590 | 394 | 3.558e-05 | 2.242e-04 | 3.478e-04 | 2.711e-02 | 5.070e-05 | 1.767e-04 | 6.464e-05 | 2.573e-03 | 4.911e-07 | 6.668e-06 | 3.087e-05 | 1.687e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 395 | 22877 | 8590 | 394 | 5.082e-05 | 3.387e-04 | 3.471e-04 | 2.709e-02 | 5.657e-05 | 1.714e-04 | 6.353e-05 | 2.583e-03 | 4.774e-07 | 6.208e-06 | 2.999e-05 | 1.689e-03 |

## III.9 Numerical controls of the transfer comparison

Steps 2..40 of the seed-1 drift sequence (three transfer events at steps 9, 19, 29) recomputed with a tighter inner loop (cap 1600 instead of 800, same events: the displacement seed is the step index) and, for SA, at half the physical step (dt 5e-4 s, 80 steps, the same three events replayed at the same physical times from the production displacement files, own half-step reference). Rows "X - production" give the change of the case-minus-reference response D_m and of the paired difference P = NN2 - GP2 between the control and the production setting at the common physical steps, as RMS over the window and relative to the production RMS (`rans_control_sensitivity.csv`, series in `rans_control_series.csv`, `RANS_Controls.png`).

| model | setting | quantity | coeff | steps | RMS | max | end | RMS of production | RMS change / production RMS |
|---|---|---|---|---|---|---|---|---|---|
| SA | production | GP2 | CL | 39 | 2.6287e-04 | 5.0585e-04 | 2.7047e-04 |  |  |
| SA | production | GP2 | CD | 39 | 6.5945e-05 | 1.0976e-04 | 7.8642e-05 |  |  |
| SA | production | GP2 | CDp | 39 | 6.5816e-05 | 1.0979e-04 | 7.8486e-05 |  |  |
| SA | production | GP2 | CDv | 39 | 2.1156e-07 | 5.7426e-07 | 1.5665e-07 |  |  |
| SA | production | NN2 | CL | 39 | 2.5978e-04 | 5.1050e-04 | 2.4323e-04 |  |  |
| SA | production | NN2 | CD | 39 | 6.6334e-05 | 1.1062e-04 | 7.5537e-05 |  |  |
| SA | production | NN2 | CDp | 39 | 6.6213e-05 | 1.1072e-04 | 7.5462e-05 |  |  |
| SA | production | NN2 | CDv | 39 | 2.1770e-07 | 5.9668e-07 | 7.6795e-08 |  |  |
| SA | production | P | CL | 39 | 1.2675e-05 | 4.6771e-05 | -2.7235e-05 |  |  |
| SA | production | P | CD | 39 | 2.2158e-06 | 6.0325e-06 | -3.1041e-06 |  |  |
| SA | production | P | CDp | 39 | 2.1982e-06 | 5.9599e-06 | -3.0242e-06 |  |  |
| SA | production | P | CDv | 39 | 3.1953e-08 | 7.9852e-08 | -7.9852e-08 |  |  |
| SA | cap1600 | GP2 | CL | 38 | 2.6244e-04 | 5.0617e-04 | 3.6019e-04 |  |  |
| SA | cap1600 | GP2 | CD | 38 | 6.5526e-05 | 1.0953e-04 | 9.5050e-05 |  |  |
| SA | cap1600 | GP2 | CDp | 38 | 6.5395e-05 | 1.0955e-04 | 9.5066e-05 |  |  |
| SA | cap1600 | GP2 | CDv | 38 | 2.1542e-07 | 5.8187e-07 | -9.7577e-09 |  |  |
| SA | cap1600 | NN2 | CL | 38 | 2.6026e-04 | 5.1050e-04 | 3.1714e-04 |  |  |
| SA | cap1600 | NN2 | CD | 38 | 6.6063e-05 | 1.1066e-04 | 8.9998e-05 |  |  |
| SA | cap1600 | NN2 | CDp | 38 | 6.5943e-05 | 1.1075e-04 | 9.0073e-05 |  |  |
| SA | cap1600 | NN2 | CDv | 38 | 2.1699e-07 | 5.8762e-07 | -6.9659e-08 |  |  |
| SA | cap1600 | P | CL | 38 | 1.1252e-05 | 4.3056e-05 | -4.3056e-05 |  |  |
| SA | cap1600 | P | CD | 38 | 2.1140e-06 | 5.1427e-06 | -5.0525e-06 |  |  |
| SA | cap1600 | P | CDp | 38 | 2.1078e-06 | 5.1134e-06 | -4.9928e-06 |  |  |
| SA | cap1600 | P | CDv | 38 | 2.6849e-08 | 7.2289e-08 | -5.9901e-08 |  |  |
| SA | finemesh | GP2 | CL | 38 | 1.0765e-04 | 2.6122e-04 | 2.6122e-04 |  |  |
| SA | finemesh | GP2 | CD | 38 | 5.3847e-05 | 8.5073e-05 | 8.5073e-05 |  |  |
| SA | finemesh | GP2 | CDp | 38 | 5.3771e-05 | 8.4468e-05 | 8.4468e-05 |  |  |
| SA | finemesh | GP2 | CDv | 38 | 2.2198e-07 | 6.0564e-07 | 6.0564e-07 |  |  |
| SA | finemesh | NN2 | CL | 38 | 1.0661e-04 | 2.6241e-04 | 2.6241e-04 |  |  |
| SA | finemesh | NN2 | CD | 38 | 5.3913e-05 | 8.6180e-05 | 8.6180e-05 |  |  |
| SA | finemesh | NN2 | CDp | 38 | 5.3834e-05 | 8.5590e-05 | 8.5590e-05 |  |  |
| SA | finemesh | NN2 | CDv | 38 | 2.1434e-07 | 5.9023e-07 | 5.9023e-07 |  |  |
| SA | finemesh | P | CL | 38 | 8.7701e-06 | 2.6501e-05 | 1.1834e-06 |  |  |
| SA | finemesh | P | CD | 38 | 9.9421e-07 | 2.9016e-06 | 1.1066e-06 |  |  |
| SA | finemesh | P | CDp | 38 | 9.9710e-07 | 2.9327e-06 | 1.1227e-06 |  |  |
| SA | finemesh | P | CDv | 38 | 1.8483e-08 | 4.4873e-08 | -1.5413e-08 |  |  |
| SA | cap1600 - production | GP2 | CL | 38 | 6.5316e-07 | 3.2314e-06 | -8.4570e-07 | 2.6267e-04 | 0.002 |
| SA | cap1600 - production | GP2 | CD | 38 | 1.3316e-07 | 7.2450e-07 | -7.2450e-07 | 6.5577e-05 | 0.002 |
| SA | cap1600 - production | GP2 | CDp | 38 | 1.3346e-07 | 7.1969e-07 | -7.1969e-07 | 6.5450e-05 | 0.002 |
| SA | cap1600 - production | GP2 | CDv | 38 | 7.6850e-09 | 1.8181e-08 | -4.7932e-09 | 2.1281e-07 | 0.036 |
| SA | cap1600 - production | NN2 | CL | 38 | 8.7936e-07 | 3.0403e-06 | 2.8687e-06 | 2.6020e-04 | 0.003 |
| SA | cap1600 - production | NN2 | CD | 38 | 6.6917e-08 | 2.5544e-07 | 2.5544e-07 | 6.6075e-05 | 0.001 |
| SA | cap1600 - production | NN2 | CDp | 38 | 6.3585e-08 | 2.4741e-07 | 2.4741e-07 | 6.5952e-05 | 0.001 |
| SA | cap1600 - production | NN2 | CDv | 38 | 5.2881e-09 | 2.0359e-08 | 8.8128e-09 | 2.2020e-07 | 0.024 |
| SA | cap1600 - production | P | CL | 38 | 1.3177e-06 | 4.6899e-06 | 3.7144e-06 | 1.2056e-05 | 0.109 |
| SA | cap1600 - production | P | CD | 38 | 1.8343e-07 | 9.7994e-07 | 9.7994e-07 | 2.1875e-06 | 0.084 |
| SA | cap1600 - production | P | CDp | 38 | 1.8186e-07 | 9.6710e-07 | 9.6710e-07 | 2.1722e-06 | 0.084 |
| SA | cap1600 - production | P | CDv | 38 | 1.1005e-08 | 3.0923e-08 | 1.3606e-08 | 2.9666e-08 | 0.371 |
| SA | finemesh - production | GP2 | CL | 38 | 1.9708e-04 | 4.5141e-04 | -9.9818e-05 | 2.6267e-04 | 0.750 |
| SA | finemesh - production | GP2 | CD | 38 | 2.1047e-05 | 4.5072e-05 | -1.0701e-05 | 6.5577e-05 | 0.321 |
| SA | finemesh - production | GP2 | CDp | 38 | 2.0945e-05 | 4.4622e-05 | -1.1318e-05 | 6.5450e-05 | 0.320 |
| SA | finemesh - production | GP2 | CDv | 38 | 3.4299e-07 | 7.5635e-07 | 6.1061e-07 | 2.1281e-07 | 1.612 |
| SA | finemesh - production | NN2 | CL | 38 | 1.9859e-04 | 4.4984e-04 | -5.1864e-05 | 2.6020e-04 | 0.763 |
| SA | finemesh - production | NN2 | CD | 38 | 2.1863e-05 | 4.8453e-05 | -3.5623e-06 | 6.6075e-05 | 0.331 |
| SA | finemesh - production | NN2 | CDp | 38 | 2.1754e-05 | 4.7993e-05 | -4.2353e-06 | 6.5952e-05 | 0.330 |
| SA | finemesh - production | NN2 | CDv | 38 | 3.6001e-07 | 7.8750e-07 | 6.6870e-07 | 2.2020e-07 | 1.635 |
| SA | finemesh - production | P | CL | 38 | 1.4281e-05 | 4.7954e-05 | 4.7954e-05 | 1.2056e-05 | 1.185 |
| SA | finemesh - production | P | CD | 38 | 2.2598e-06 | 7.1391e-06 | 7.1391e-06 | 2.1875e-06 | 1.033 |
| SA | finemesh - production | P | CDp | 38 | 2.2554e-06 | 7.0826e-06 | 7.0826e-06 | 2.1722e-06 | 1.038 |
| SA | finemesh - production | P | CDv | 38 | 2.6429e-08 | 7.1203e-08 | 5.8095e-08 | 2.9666e-08 | 0.891 |
| SST | production | GP2 | CL | 39 | 1.6331e-04 | 3.9496e-04 | -3.6097e-05 |  |  |
| SST | production | GP2 | CD | 39 | 5.4299e-05 | 8.7759e-05 | 4.6261e-05 |  |  |
| SST | production | GP2 | CDp | 39 | 5.5066e-05 | 8.9198e-05 | 4.7708e-05 |  |  |
| SST | production | GP2 | CDv | 39 | 8.3380e-07 | 1.4587e-06 | -1.4463e-06 |  |  |
| SST | production | NN2 | CL | 39 | 1.6329e-04 | 3.8733e-04 | -1.6240e-05 |  |  |
| SST | production | NN2 | CD | 39 | 5.4934e-05 | 8.9994e-05 | 4.7353e-05 |  |  |
| SST | production | NN2 | CDp | 39 | 5.5692e-05 | 9.1425e-05 | 4.8838e-05 |  |  |
| SST | production | NN2 | CDv | 39 | 8.2698e-07 | 1.4849e-06 | -1.4849e-06 |  |  |
| SST | production | P | CL | 39 | 6.7012e-06 | 1.9857e-05 | 1.9857e-05 |  |  |
| SST | production | P | CD | 39 | 1.5769e-06 | 3.6014e-06 | 1.0918e-06 |  |  |
| SST | production | P | CDp | 39 | 1.5744e-06 | 3.5803e-06 | 1.1299e-06 |  |  |
| SST | production | P | CDv | 39 | 1.6001e-08 | 3.8557e-08 | -3.8557e-08 |  |  |
| SST | cap1600 | GP2 | CL | 38 | 1.6485e-04 | 3.8718e-04 | 5.5082e-05 |  |  |
| SST | cap1600 | GP2 | CD | 38 | 5.4579e-05 | 8.8809e-05 | 6.2680e-05 |  |  |
| SST | cap1600 | GP2 | CDp | 38 | 5.5478e-05 | 9.0540e-05 | 6.4234e-05 |  |  |
| SST | cap1600 | GP2 | CDv | 38 | 9.5799e-07 | 1.7269e-06 | -1.5514e-06 |  |  |
| SST | cap1600 | NN2 | CL | 38 | 1.6540e-04 | 3.8150e-04 | 5.3618e-05 |  |  |
| SST | cap1600 | NN2 | CD | 38 | 5.5229e-05 | 9.1087e-05 | 6.2711e-05 |  |  |
| SST | cap1600 | NN2 | CDp | 38 | 5.6119e-05 | 9.2832e-05 | 6.4309e-05 |  |  |
| SST | cap1600 | NN2 | CDv | 38 | 9.5397e-07 | 1.7475e-06 | -1.5964e-06 |  |  |
| SST | cap1600 | P | CL | 38 | 5.6821e-06 | 1.6198e-05 | -1.4634e-06 |  |  |
| SST | cap1600 | P | CD | 38 | 1.6213e-06 | 3.6526e-06 | 3.0618e-08 |  |  |
| SST | cap1600 | P | CDp | 38 | 1.6187e-06 | 3.6230e-06 | 7.5436e-08 |  |  |
| SST | cap1600 | P | CDv | 38 | 1.8193e-08 | 4.4952e-08 | -4.4952e-08 |  |  |
| SST | cap6400 | GP2 | CL | 38 | 1.6633e-04 | 3.8652e-04 | 4.4478e-05 |  |  |
| SST | cap6400 | GP2 | CD | 38 | 5.4743e-05 | 9.0011e-05 | 6.1807e-05 |  |  |
| SST | cap6400 | GP2 | CDp | 38 | 5.5733e-05 | 9.1911e-05 | 6.3250e-05 |  |  |
| SST | cap6400 | GP2 | CDv | 38 | 1.0507e-06 | 1.9288e-06 | -1.4411e-06 |  |  |
| SST | cap6400 | NN2 | CL | 38 | 1.6424e-04 | 3.6871e-04 | 3.8301e-05 |  |  |
| SST | cap6400 | NN2 | CD | 38 | 5.5271e-05 | 9.2059e-05 | 6.1404e-05 |  |  |
| SST | cap6400 | NN2 | CDp | 38 | 5.6247e-05 | 9.3990e-05 | 6.2895e-05 |  |  |
| SST | cap6400 | NN2 | CDv | 38 | 1.0435e-06 | 1.9495e-06 | -1.4910e-06 |  |  |
| SST | cap6400 | P | CL | 38 | 8.3400e-06 | 2.7610e-05 | -6.1769e-06 |  |  |
| SST | cap6400 | P | CD | 38 | 1.6366e-06 | 3.9253e-06 | -4.0340e-07 |  |  |
| SST | cap6400 | P | CDp | 38 | 1.6303e-06 | 3.8783e-06 | -3.5440e-07 |  |  |
| SST | cap6400 | P | CDv | 38 | 3.3980e-08 | 8.7479e-08 | -4.9893e-08 |  |  |
| SST | cap1600 - production | GP2 | CL | 38 | 6.1744e-06 | 1.3551e-05 | -8.0065e-06 | 1.6534e-04 | 0.037 |
| SST | cap1600 - production | GP2 | CD | 38 | 6.1471e-07 | 1.5461e-06 | -7.3941e-07 | 5.4494e-05 | 0.011 |
| SST | cap1600 - production | GP2 | CDp | 38 | 7.2612e-07 | 1.8351e-06 | -6.4663e-07 | 5.5247e-05 | 0.013 |
| SST | cap1600 - production | GP2 | CDv | 38 | 2.0382e-07 | 3.7229e-07 | -9.2768e-08 | 8.1146e-07 | 0.251 |
| SST | cap1600 - production | NN2 | CL | 38 | 6.0004e-06 | 1.3805e-05 | -6.9360e-06 | 1.6540e-04 | 0.036 |
| SST | cap1600 - production | NN2 | CD | 38 | 5.7343e-07 | 1.3357e-06 | -6.8166e-07 | 5.5120e-05 | 0.010 |
| SST | cap1600 - production | NN2 | CDp | 38 | 6.9920e-07 | 1.6420e-06 | -5.6971e-07 | 5.5861e-05 | 0.013 |
| SST | cap1600 - production | NN2 | CDv | 38 | 2.0719e-07 | 3.7987e-07 | -1.1241e-07 | 8.0241e-07 | 0.258 |
| SST | cap1600 - production | P | CL | 38 | 8.1872e-07 | 2.2582e-06 | 1.0705e-06 | 5.9760e-06 | 0.137 |
| SST | cap1600 - production | P | CD | 38 | 9.5874e-08 | 2.6040e-07 | 5.7745e-08 | 1.5877e-06 | 0.060 |
| SST | cap1600 - production | P | CDp | 38 | 9.6319e-08 | 2.8537e-07 | 7.6920e-08 | 1.5845e-06 | 0.061 |
| SST | cap1600 - production | P | CDv | 38 | 9.4598e-09 | 2.4150e-08 | -1.9645e-08 | 1.4955e-08 | 0.633 |
| SST | cap6400 - production | GP2 | CL | 38 | 1.1461e-05 | 2.8733e-05 | -1.8610e-05 | 1.6534e-04 | 0.069 |
| SST | cap6400 - production | GP2 | CD | 38 | 1.1599e-06 | 3.1448e-06 | -1.6124e-06 | 5.4494e-05 | 0.021 |
| SST | cap6400 - production | GP2 | CDp | 38 | 1.3620e-06 | 3.6619e-06 | -1.6306e-06 | 5.5247e-05 | 0.025 |
| SST | cap6400 - production | GP2 | CDv | 38 | 3.2902e-07 | 6.2197e-07 | 1.7592e-08 | 8.1146e-07 | 0.405 |
| SST | cap6400 - production | NN2 | CL | 38 | 1.4158e-05 | 3.1522e-05 | -2.2253e-05 | 1.6540e-04 | 0.086 |
| SST | cap6400 - production | NN2 | CD | 38 | 1.3957e-06 | 3.3941e-06 | -1.9886e-06 | 5.5120e-05 | 0.025 |
| SST | cap6400 - production | NN2 | CDp | 38 | 1.6020e-06 | 3.9756e-06 | -1.9835e-06 | 5.5861e-05 | 0.029 |
| SST | cap6400 - production | NN2 | CDv | 38 | 3.3536e-07 | 6.2751e-07 | -6.9941e-09 | 8.0241e-07 | 0.418 |
| SST | cap6400 - production | P | CL | 38 | 4.0682e-06 | 1.0735e-05 | -3.6430e-06 | 5.9760e-06 | 0.681 |
| SST | cap6400 - production | P | CD | 38 | 4.0737e-07 | 1.0936e-06 | -3.7627e-07 | 1.5877e-06 | 0.257 |
| SST | cap6400 - production | P | CDp | 38 | 4.2113e-07 | 1.1646e-06 | -3.5291e-07 | 1.5845e-06 | 0.266 |
| SST | cap6400 - production | P | CDv | 38 | 2.6916e-08 | 7.5739e-08 | -2.4586e-08 | 1.4955e-08 | 1.800 |

Inner-loop record of the control runs (`rans_control_inner.csv`):

| model | setting | method | steps | cap | steps at cap | median inner | max inner | median log10 rms(rho) at end | worst |
|---|---|---|---|---|---|---|---|---|---|
| SA | production | GP2 | 398 | 800 | 144 | 697 | 799 | -12.00 | -9.12 |
| SA | production | NN2 | 398 | 800 | 156 | 711 | 799 | -12.00 | -9.21 |
| SA | cap1600 | GP2 | 38 | 1600 | 0 | 575 | 957 | -12.00 | -12.00 |
| SA | cap1600 | NN2 | 38 | 1600 | 0 | 593 | 1044 | -12.00 | -12.00 |
| SA | finemesh | GP2 | 38 | 800 | 7 | 483 | 799 | -12.00 | -10.95 |
| SA | finemesh | NN2 | 38 | 800 | 5 | 466 | 799 | -12.00 | -10.90 |
| SST | production | GP2 | 398 | 800 | 333 | 799 | 799 | -11.10 | -10.11 |
| SST | production | NN2 | 398 | 800 | 323 | 799 | 799 | -11.24 | -9.52 |
| SST | cap1600 | GP2 | 38 | 1600 | 28 | 1599 | 1599 | -11.33 | -10.35 |
| SST | cap1600 | NN2 | 38 | 1600 | 28 | 1599 | 1599 | -11.33 | -10.46 |
| SST | cap6400 | GP2 | 38 | 6400 | 28 | 6399 | 6399 | -11.38 | -10.29 |
| SST | cap6400 | NN2 | 38 | 6400 | 30 | 6399 | 6399 | -11.28 | -10.38 |

Event replay check (`rans_control_replay.csv`): 30 displacement files compared with the production case by md5, 30 identical, 0 different.

## III.10 Viscous-mesh sensitivity of the no-transfer baseline

`Meshes/Coarse_Sharp_Finer.su2` (32102 nodes, 15128 quadrilaterals, 312 airfoil nodes, same first off-wall spacing 5e-6, finer wall-normal growth and tangential spacing) against the production grid `Meshes/Coarse_Sharp.su2` (23271 nodes): steady state and 800-step no-transfer reference per closure, no transfer involved (`rans_fine_mesh.csv`, `rans_fine_surface.csv`, `RANS_FineMesh_<MODEL>.png`).

| model | quantity | coeff | coarse | fine | fine - coarse | % |
|---|---|---|---|---|---|---|
| SA | steady | CL | 0.676257 | 0.670838 | -5.419e-03 | -0.801 |
| SA | steady | CD | 0.010229 | 0.010048 | -1.812e-04 | -1.771 |
| SA | steady | CDp | 0.003559 | 0.003401 | -1.580e-04 | -4.439 |
| SA | steady | CDv | 0.006671 | 0.006647 | -2.400e-05 | -0.360 |
| SA | steady | CMz | 0.000671 | -0.000530 | -1.201e-03 | -179.066 |
| SA | reference mean 100..799 | CL | 0.676257 | 0.670838 | -5.419e-03 | -0.801 |
| SA | reference range 100..799 | CL | 6.40e-09 | 1.12e-07 |  |  |
| SA | reference mean 100..799 | CD | 0.010229 | 0.010048 | -1.812e-04 | -1.771 |
| SA | reference range 100..799 | CD | 3.94e-09 | 2.38e-09 |  |  |
| SA | fine reference inner iterations min..max |  | 10 | 10 |  |  |
| SST | steady | CL | 0.675605 | 0.668596 | -7.008e-03 | -1.037 |
| SST | steady | CD | 0.009858 | 0.009705 | -1.531e-04 | -1.553 |
| SST | steady | CDp | 0.003469 | 0.003321 | -1.480e-04 | -4.266 |
| SST | steady | CDv | 0.006389 | 0.006385 | -4.000e-06 | -0.063 |
| SST | steady | CMz | 0.000382 | -0.001155 | -1.537e-03 | -402.191 |
| SST | reference mean 100..799 | CL | 0.675603 | 0.668596 | -7.007e-03 | -1.037 |
| SST | reference range 100..799 | CL | 7.66e-07 | 3.81e-08 |  |  |
| SST | reference mean 100..799 | CD | 0.009858 | 0.009705 | -1.531e-04 | -1.553 |
| SST | reference range 100..799 | CD | 5.89e-08 | 2.09e-09 |  |  |
| SST | fine reference inner iterations min..max |  | 56 | 80 |  |  |

Three-level sequence (`rans_mesh_levels.csv`; `Meshes/Coarse_Sharp_FinerFiner.su2`: 48228 nodes, 23764 quadrilaterals, 516 airfoil nodes, same first spacing): steady loads per level with the effective spacing ratio sqrt(N_coarse/N), and the observed order from the three levels where the sequence is monotone (Richardson extrapolation of the finest value):

| model | level | nodes | h / h_coarse | CL | CD | CDp | CDv |
|---|---|---|---|---|---|---|---|
| SA | coarse | 23271 | 1.0000 | 0.676257 | 0.010229 | 0.003559 | 0.006671 |
| SA | fine | 32102 | 0.8514 | 0.670838 | 0.010048 | 0.003401 | 0.006647 |
| SA | finer | 48228 | 0.6946 | 0.668535 | 0.009996 | 0.003330 | 0.006665 |
| SA | observed order CL |  |  | 4.21 | extrapolated 0.666833 | coarse-to-finer change 7.722e-03 | finer error estimate 1.701e-03 |
| SA | observed order CD |  |  | 6.07 | extrapolated 0.009974 | coarse-to-finer change 2.338e-04 | finer error estimate 2.154e-05 |
| SA | observed order CDp |  |  | 3.93 | extrapolated 0.003272 | coarse-to-finer change 2.290e-04 | finer error estimate 5.794e-05 |
| SA | observed order CDv |  |  | not monotone or ratios differ | 2.400e-05 | -1.800e-05 |  |
| SST | coarse | 23271 | 1.0000 | 0.675605 | 0.009858 | 0.003469 | 0.006389 |
| SST | fine | 32102 | 0.8514 | 0.668596 | 0.009705 | 0.003321 | 0.006385 |
| SST | finer | 48228 | 0.6946 | 0.666394 | 0.009644 | 0.003247 | 0.006397 |
| SST | observed order CL |  |  | 5.69 | extrapolated 0.665385 | coarse-to-finer change 9.211e-03 | finer error estimate 1.009e-03 |
| SST | observed order CD |  |  | 4.49 | extrapolated 0.009603 | coarse-to-finer change 2.144e-04 | finer error estimate 4.101e-05 |
| SST | observed order CDp |  |  | 3.41 | extrapolated 0.003173 | coarse-to-finer change 2.220e-04 | finer error estimate 7.400e-05 |
| SST | observed order CDv |  |  | not monotone or ratios differ | 4.000e-06 | -1.200e-05 |  |

Surface pressure and skin friction, fine field interpolated onto the coarse airfoil nodes along each side:

| model | side | coarse nodes | fine nodes | RMS dCp | max dCp | x | RMS dCf | max dCf | x | max y+ coarse | max y+ fine |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SA | upper | 129 | 156 | 1.044e-02 | 4.909e-02 | 1.000 | 1.448e-04 | 9.646e-04 | 0.998 | 2.765 | 2.750 |
| SA | lower | 127 | 156 | 2.752e-03 | 8.185e-03 | 0.000 | 9.767e-05 | 8.015e-04 | 0.998 | 2.406 | 2.390 |
| SST | upper | 129 | 156 | 1.246e-02 | 4.457e-02 | 1.000 | 1.431e-04 | 9.512e-04 | 0.998 | 2.698 | 2.691 |
| SST | lower | 127 | 156 | 4.081e-03 | 1.244e-02 | 0.000 | 9.341e-05 | 8.005e-04 | 0.998 | 2.281 | 2.266 |

## III.7b Cost of one transfer event

Wall-clock times from the launcher's per-event timing files (`rans_timing.csv`). "Transfer" is the sum of the two interpolations of an event (state k and BDF2 level k−1), including geometry search and supermesh construction on each call (nothing is reused between events); the first event includes the numba just-in-time compilation. Mesh move = RandomMeshMovements. CFD = one 10-step window on 4 MPI ranks. Transfer and mesh motion ran on 4 threads. No memory measurement was recorded.

| case | model | method | transfer, first event (s) | transfer per event, median (s) | min | max | mesh move per event (s) | CFD per 10-step window (s) |
|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | SA | ConsGalerkinProj | 3.38 | 7.33 | 6.58 | 13.56 | 2.01 | 19.41 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | SA | ConsGalerkinProj | 2.52 | 5.07 | 4.90 | 5.79 | 1.47 | 10.46 |
| SA 0.01-20 NN Drift | SA | NN | 2.80 | 3.54 | 3.48 | 4.74 | 1.45 | 10.11 |
| SA 0.01-20 NN Drift Seed2 | SA | NN | 1.80 | 4.87 | 3.56 | 5.71 | 1.94 | 19.16 |
| SST 0.01-20 ConsGalerkinProj Drift | SST | ConsGalerkinProj | 2.73 | 5.34 | 5.00 | 6.56 | 1.51 | 11.68 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | SST | ConsGalerkinProj | 3.23 | 6.11 | 5.49 | 8.10 | 1.67 | 13.19 |
| SST 0.01-20 NN Drift | SST | NN | 2.21 | 5.52 | 4.14 | 7.54 | 2.16 | 21.27 |
| SST 0.01-20 NN Drift Seed2 | SST | NN | 2.02 | 4.43 | 4.03 | 6.21 | 1.78 | 14.49 |

## III.7 Figures and files

Figures (`Figures/RANS/`):

- `RANS_Controls.png`
- `RANS_FineMesh_SA.png`
- `RANS_FineMesh_SST.png`
- `RANS_SA_AcousticDelay_CL.png`
- `RANS_SA_DeformedReference.png`
- `RANS_SA_Drift_DragSplit_Surface.png`
- `RANS_SA_ForceDeviation.png`
- `RANS_SA_OneShot_DragSplit.png`
- `RANS_SA_OneShot_DragSplit_Surface.png`
- `RANS_SA_PressureProbe_Pressure.png`
- `RANS_SA_SurfaceNorms.png`
- `RANS_SA_TransferIntegrals.png`
- `RANS_SA_TurbState.png`
- `RANS_SA_WD_0.01-20_ConsGalerkinProj_Drift_Second_Cf_Cp.png`
- `RANS_SA_WD_0.01-20_ConsGalerkinProj_Drift_Seed2_Second_Cf_Cp.png`
- `RANS_SA_WD_0.01-20_NN_Drift_Second_Cf_Cp.png`
- `RANS_SA_WD_0.01-20_NN_Drift_Seed2_Second_Cf_Cp.png`
- `RANS_SA_WD_0.05-0.10_ConsGalerkinProj_OneShot_Second_Cf_Cp.png`
- `RANS_SA_WD_0.10-0.15_ConsGalerkinProj_OneShot_SF2_Second_Cf_Cp.png`
- `RANS_SA_WD_0.10-0.15_ConsGalerkinProj_OneShot_Second_Cf_Cp.png`
- `RANS_SA_WD_0.20-0.25_ConsGalerkinProj_OneShot_Long_Second_Cf_Cp.png`
- `RANS_SA_WD_0.20-0.25_ConsGalerkinProj_OneShot_Second_Cf_Cp.png`
- `RANS_SST_AcousticDelay_CL.png`
- `RANS_SST_DeformedReference.png`
- `RANS_SST_Drift_DragSplit_Surface.png`
- `RANS_SST_ForceDeviation.png`
- `RANS_SST_OneShot_DragSplit.png`
- `RANS_SST_OneShot_DragSplit_Surface.png`
- `RANS_SST_PressureProbe_Pressure.png`
- `RANS_SST_SurfaceNorms.png`
- `RANS_SST_TransferIntegrals.png`
- `RANS_SST_TurbState.png`
- `RANS_SST_WD_0.01-20_ConsGalerkinProj_Drift_Second_Cf_Cp.png`
- `RANS_SST_WD_0.01-20_ConsGalerkinProj_Drift_Seed2_Second_Cf_Cp.png`
- `RANS_SST_WD_0.01-20_NN_Drift_Second_Cf_Cp.png`
- `RANS_SST_WD_0.01-20_NN_Drift_Seed2_Second_Cf_Cp.png`
- `RANS_SST_WD_0.05-0.10_ConsGalerkinProj_OneShot_Second_Cf_Cp.png`
- `RANS_SST_WD_0.10-0.15_ConsGalerkinProj_OneShot_SF2_Second_Cf_Cp.png`
- `RANS_SST_WD_0.10-0.15_ConsGalerkinProj_OneShot_Second_Cf_Cp.png`
- `RANS_SST_WD_0.20-0.25_ConsGalerkinProj_OneShot_Long_Second_Cf_Cp.png`
- `RANS_SST_WD_0.20-0.25_ConsGalerkinProj_OneShot_Second_Cf_Cp.png`
- `RANS_T04_Sensitivity.png`

Tables (`RANS_PostProcessing/`):

- `rans_code_provenance.csv`
- `rans_control_inner.csv`
- `rans_control_replay.csv`
- `rans_control_sensitivity.csv`
- `rans_control_series.csv`
- `rans_deformed_reference.csv`
- `rans_deformed_surface.csv`
- `rans_drift_series.csv`
- `rans_drift_trends.csv`
- `rans_fine_mesh.csv`
- `rans_fine_surface.csv`
- `rans_inner_convergence.csv`
- `rans_mesh_levels.csv`
- `rans_mesh_quality.csv`
- `rans_mesh_replay.csv`
- `rans_oneshot_fit.csv`
- `rans_oneshot_lags.csv`
- `rans_oneshot_series.csv`
- `rans_pressure_probe.csv`
- `rans_protected_region.csv`
- `rans_protected_summary.csv`
- `rans_runtime_units.csv`
- `rans_surface_forces.csv`
- `rans_surface_norms.csv`
- `rans_surface_profiles.csv`
- `rans_timing.csv`
- `rans_transfer_integrals.csv`
- `rans_turb_state_norms_SA.csv`
- `rans_turb_state_norms_SST.csv`

Scripts: `RunAll_RANS.sh` (cases), `Unsteady_Ref_RANS_<MODEL>_Long/run.sh` (references), `RANS_ForceAnalysis.py` (III.3, III.4), `RANS_SurfaceAnalysis.py` (III.5 and the full-precision drag split), `RANS_TurbulenceState.py` (III.6), `RANS_PressureProbe.py` (III.4b), `RANS_DeformedReference.py` (III.2b), `RANS_WriteReportData.py` (this section); `run_rans_postprocessing.sh` runs them in order.

## III.8 Notes for the writer

- Every number above is read from the run outputs; the tables are generated, not typed.
- `forces_breakdown` carries 6 decimals; CDp/CDv from it are at rounding when the deviation is below ~2e-6. The surface-integrated split (`rans_surface_forces.csv`) resolves ~1e-9 but uses a trapezoid quadrature, so its absolute values differ from SU2's at the 1e-4 level; only differences case − reference at the same step should be quoted from it.
- The transported-variable integrals change between transfers because of production, destruction and diffusion; only the source-to-target defect at a transfer is a property of the transfer.
- Volume and surface VTU fields are single precision (pressure ULP 6e-8, i.e. 1.2e-5 of the dynamic pressure); a force deviation of 1e-5 in CL corresponds to surface pressure changes below one ULP, so field-based probes cannot resolve the smallest force responses. History CSVs and forces_breakdown are double/6-digit.
- The pre-event floor of each one-shot case is the bit-identity check between the case and the reference; a non-zero floor means the two runs differ before any transfer and the reference must be questioned.
- Positivity of the transferred turbulence variables: see the clipping table in III.6 for which cases ran with clipFieldsToMinimum on and how many nodes it touched. The SA cases ran with it off; their largest undershoot of nu_tilde was of order 1e-16.
- Wall distance, y+ and the eddy viscosity are recomputed by the solver on every mesh; the node-wise comparison uses only nodes that never moved.

<!-- RANS-DATA-END -->

---

## 15. Scope

One geometry, one Mach number, one mesh resolution class. The drift result uses two seeds
per method, enough to separate a reproducible drag trend from a non-reproducible lift
trend but not to characterise the distribution of either. The propagation speed is
measured on a single event; the d = 0.10 case corroborates it but has only three shells
outside the near-field exclusion, so its fit is a consistency check rather than an
independent determination.

For Part II: one mesh pair and one transfer direction. The pair shares its boundary
discretisation, so it does not cover mesh pairs whose surfaces are re-discretised under
surface adaptation. The flow is inviscid, so no turbulence or transition variables are
transferred and no skin friction is reported. The steady experiment applies the same
spatial defect at both history levels; a general pair of non-collinear history defects
requires a weakly unsteady phase sweep. The linearised prediction from the target
Jacobian is not evaluated, as the solver exposes no ready route to it.

## 16. Reproduction

Part I:

    RunAll.sh                     builds and runs the drift and one-shot cases
    ForceDriftAnalysis.py         section 4
    AcousticDelayAnalysis.py      section 5
    PressureProbeAnalysis.py      section 6

Part II, under `BDF2/`:

    RunTransfers.sh               transfers U_A* onto mesh B for all 12 methods
    RunHistories.py               builds and runs the history matrix and residual probes
    bdf2_analysis.py              produces the CSV deliverables
    bdf2_utils.py                 restart I/O, nodal volumes, volume-weighted norms
    GP3_diagnosis/                the patch, coverage and clipper tests behind section 14

Case directories are skipped when already complete, so any stage can be resumed; pass
`FORCE=1` to rebuild. Figures are written to `Figures/`: `ForceHistories.png` and
`ForceDeviation.png` for the drift, `AcousticDelay_CL.png` for the force-based delay, and
`PressureProbe_Pressure.png` for the shell arrivals and the speed fit. Tabular outputs are
`naca_native_convergence.csv`, `naca_transfer_defects.csv`,
`naca_initial_residual_identity.csv`, `naca_first_step_response.csv`,
`naca_superposition.csv`, `naca_force_history.csv` and `naca_surface_delta_cp.csv`.
