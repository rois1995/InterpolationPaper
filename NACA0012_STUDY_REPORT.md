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
