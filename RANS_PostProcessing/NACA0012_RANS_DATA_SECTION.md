<!-- RANS-DATA-BEGIN -->

# Part III — RANS extension: data for the results section

Same experiment as Part I on a wall-resolved hybrid grid with the Spalart–Allmaras (SA) and Menter SST turbulence models. Two additions against the inviscid study: the transported turbulence state is part of what is transferred, and the drag is reported split into its pressure and friction parts. This part collects the measured quantities; it does not interpret them.

## III.1 Setup

**Grid.** `Meshes/Coarse.su2`: 29956 nodes, 14070 quadrilaterals in the wall layer and 31339 triangles outside it; 370 airfoil nodes (370 surface elements), farfield at 100 chords. First off-wall spacing 5.00e-06 to 5.14e-06; the quadrilateral layer extends to 0.053 chords from the wall. The trailing edge is blunt (five nodes at x = 1). The same grid, undeformed, is used by both models and by the references.

**SA solver settings** (`Base_Interp_RANS_SA/Unsteady_SecondOrder.cfg`): SOLVER = RANS, KIND_TURB_MODEL = SA, TIME_MARCHING = DUAL_TIME_STEPPING-2ND_ORDER, TIME_STEP = 1e-3, INNER_ITER = 50, MACH_NUMBER = 0.1, AOA = 6.0, REYNOLDS_NUMBER = 6.0E6, REF_DIMENSIONALIZATION = FREESTREAM_VEL_EQ_MACH, NUM_METHOD_GRAD = WEIGHTED_LEAST_SQUARES, CFL_NUMBER = 6.0, LINEAR_SOLVER = FGMRES, LINEAR_SOLVER_PREC = ILU, LINEAR_SOLVER_ITER = 60, CONV_NUM_METHOD_FLOW = HLLC, MUSCL_FLOW = YES, SLOPE_LIMITER_FLOW = VENKATAKRISHNAN, CONV_NUM_METHOD_TURB = SCALAR_UPWIND, MUSCL_TURB = NO.

**SST solver settings** (`Base_Interp_RANS_SST/Unsteady_SecondOrder.cfg`): SOLVER = RANS, KIND_TURB_MODEL = SST, TIME_MARCHING = DUAL_TIME_STEPPING-2ND_ORDER, TIME_STEP = 1e-3, INNER_ITER = 50, MACH_NUMBER = 0.1, AOA = 6.0, REYNOLDS_NUMBER = 6.0E6, REF_DIMENSIONALIZATION = FREESTREAM_VEL_EQ_MACH, NUM_METHOD_GRAD = WEIGHTED_LEAST_SQUARES, CFL_NUMBER = 6.0, LINEAR_SOLVER = FGMRES, LINEAR_SOLVER_PREC = ILU, LINEAR_SOLVER_ITER = 60, CONV_NUM_METHOD_FLOW = HLLC, MUSCL_FLOW = YES, SLOPE_LIMITER_FLOW = VENKATAKRISHNAN, CONV_NUM_METHOD_TURB = SCALAR_UPWIND, MUSCL_TURB = NO.

**Transfer and mesh motion** (both models): InterpMethod = ConsGalerkinProj, InterpOrder = Second, WLSDegree = 1, BoundaryTreatment = False, BoundaryPointValue = False, EnforceConservation = False, ConservationCorrection = True, GalerkinLimit = True, findOneToOneCorrespondence = True, MinWallDistance = 0.01, MaxWallDistance = 20, FreezeBoundary = True, Field = "white", SafeFactor = 4.0, QualityFloor = 0.2, Seed = 12345. The launcher sets the drift band and one-shot bands per case (see the case names). The wall band `MinWallDistance` and every boundary node are held fixed; the wall-layer quadrilaterals therefore never move. The transferred fields are Density, Momentum_x, Momentum_y, Energy and the model's transported turbulence variables (Nu_Tilde for SA; Turb_Kin_Energy and Omega for SST, both transferred as rho-weighted quantities).

**Non-dimensionalisation.** `REF_DIMENSIONALIZATION = FREESTREAM_VEL_EQ_MACH`: the free-stream sound speed is the velocity reference, so a = 1 and V_inf = 0.1 in code units. With dt = 1e-3 an acoustic front travels 1e-3 chords per step and the predicted force-response lag for a perturbation launched at radius d is **1000·d** steps (Euler study: 845·d).

**Reference runs.** `Unsteady_Ref_RANS_<MODEL>_Long`: 800 steps on the undeformed grid with no transfer, started from the same two restart levels as every case (the converged steady state and the first-order start-up state written by the launcher). Cases are compared to it step by step.

## III.2 Steady solutions and references

Converged steady state (`Steady_RANS_<MODEL>/forces_breakdown.dat`, limiter on, 30 000 iterations, residual of density below 1e-8):

| model | CL | CD | CD pressure | CD friction | CMz |
|---|---|---|---|---|---|
| SA | 0.655228 | 0.008428 | 0.001786 | 0.006642 | -0.003895 |
| SST | 0.653238 | 0.008120 | 0.001751 | 0.006369 | -0.004481 |

Reference (no transfer) force histories:

| model | steps | CL first | CL step 399 | CL last | CD first | CD last | CL std | CD std |
|---|---|---|---|---|---|---|---|---|
| SA | 2..799 | 0.655218 | 0.654948 | 0.654927 | 0.008424 | 0.008423 | 6.397e-05 | 5.637e-06 |
| SST | 2..799 | 0.653232 | 0.653002 | 0.652943 | 0.008122 | 0.008123 | 6.911e-05 | 4.950e-06 |

## III.2b Reference on a moved grid

The first randomly moved grid of the drift cases (`Meshes/Coarse_moved_00009.su2`, the same for both models: 15681 of 29956 nodes displaced, wall layer and boundary fixed) is run without any transfer: the steady RANS converged on it from the original-grid steady state, then 400 unsteady steps from that state. Differences to the same runs on the original grid are the effect of the grid alone. The drift cases compute their first window, steps 10..19, on exactly this grid from the transferred state; their deviation in that window is listed for comparison (`rans_deformed_reference.csv`, `RANS_<MODEL>_DeformedReference.png`).

| model | what | n | dCL | dCD | dCDp | dCDv | dCMz |
|---|---|---|---|---|---|---|---|
| SA | steady | 14999 | -1.891e-03 | -1.189e-05 | -1.000e-05 | -3.000e-06 | -3.284e-04 |
| SA | unsteady 100..399 | 300 | -1.941e-03 (trend -1.7e-05 ± 4.8e-06) | -1.425e-05 (trend -2.2e-07 ± 5.8e-07) |  |  | -3.736e-04 (trend -9.9e-06 ± 2.2e-06) |
| SA | drift case first window 10..19: RANS_SA_WD_0.01-20_ConsGalerkinProj_Drift_Second | 10 | 1.318e-05 (min -2.3e-05, max 3.9e-05) | 7.538e-06 (min -4.7e-06, max 2.4e-05) |  |  | -1.541e-05 (min -2.5e-05, max 8.2e-06) |
| SA | drift case first window 10..19: RANS_SA_WD_0.01-20_NN_Drift_Second | 10 | 1.550e-05 (min -2.6e-05, max 4.6e-05) | 6.469e-06 (min -4.4e-06, max 2.0e-05) |  |  | -1.510e-05 (min -2.5e-05, max 8.4e-06) |
| SA | steady, surface integration | 0 | -1.891e-03 | -1.188e-05 | -9.427e-06 | -2.457e-06 |  |
| SST | steady | 14999 | -9.649e-04 | -9.910e-07 | -1.000e-06 | 1.000e-06 | -1.422e-04 |
| SST | unsteady 100..399 | 300 | -1.022e-03 (trend -1.2e-05 ± 5.6e-06) | -1.420e-05 (trend -9.5e-08 ± 5.1e-07) |  |  | -2.057e-04 (trend -1.0e-05 ± 2.5e-06) |
| SST | drift case first window 10..19: RANS_SST_WD_0.01-20_ConsGalerkinProj_Drift_Second | 10 | 1.150e-05 (min -2.6e-05, max 4.1e-05) | 7.748e-06 (min -4.7e-06, max 2.5e-05) |  |  | -1.723e-05 (min -3.6e-05, max 9.7e-06) |
| SST | drift case first window 10..19: RANS_SST_WD_0.01-20_NN_Drift_Second | 10 | 1.406e-05 (min -2.9e-05, max 4.8e-05) | 6.682e-06 (min -4.4e-06, max 2.0e-05) |  |  | -1.675e-05 (min -3.4e-05, max 1.0e-05) |
| SST | steady, surface integration | 0 | -9.649e-04 | -1.004e-06 | -1.396e-06 | 3.922e-07 |  |

Surface difference between the two steady states, node to node on the 370 airfoil nodes:

| model | RMS dCp | max dCp | x of max | RMS dCf_x | max dCf_x | x of max |
|---|---|---|---|---|---|---|
| SA | 6.126e-03 | 2.919e-02 | 0.3821 | 2.467e-05 | 1.232e-04 | 0.1019 |
| SST | 6.351e-03 | 2.595e-02 | 0.1773 | 2.145e-05 | 1.115e-04 | 0.1019 |

## III.3 Force drift under repeated transfer (transfer every 10 steps, 400 steps)

Deviation from the reference, OLS trend per 100 steps with the standard error of the slope. CL, CD, CMz from the history files (full precision); CDp, CDv from `forces_breakdown` (6 decimals; "rounding" marks a signal below that precision). "resolved" = |trend| > 2 standard errors.

| case | coeff | ref mean | offset | scatter | max |dev| | trend /100it | ± stderr | resolved | % /100it |
|---|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | CL | 6.550247e-01 | -4.926e-05 | 1.393e-04 | 4.482e-04 | -4.750e-05 | 5.61e-06 | yes | -0.0073 |
| SA 0.01-20 ConsGalerkinProj Drift | CD | 8.427870e-03 | 1.042e-05 | 3.098e-05 | 1.011e-04 | -1.750e-06 | 1.35e-06 | no | -0.0208 |
| SA 0.01-20 ConsGalerkinProj Drift | CDp | 1.785618e-03 | 1.081e-05 | 3.094e-05 | 1.020e-04 | -1.350e-06 | 1.35e-06 | no | -0.0756 |
| SA 0.01-20 ConsGalerkinProj Drift | CDv | 6.642000e-03 | -2.739e-07 | 4.459e-07 | 1.000e-06 | -2.996e-07 | 1.24e-08 | yes (rounding) | -0.0045 |
| SA 0.01-20 ConsGalerkinProj Drift | CMz | -3.954594e-03 | -2.431e-05 | 6.297e-05 | 1.916e-04 | -1.613e-05 | 2.63e-06 | yes | 0.4078 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | CL | 6.550247e-01 | -5.589e-05 | 1.467e-04 | 4.038e-04 | -5.840e-05 | 5.71e-06 | yes | -0.0089 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | CD | 8.427870e-03 | 6.115e-06 | 2.943e-05 | 1.123e-04 | -1.704e-06 | 1.28e-06 | no | -0.0202 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | CDp | 1.785618e-03 | 6.606e-06 | 2.942e-05 | 1.120e-04 | -1.248e-06 | 1.29e-06 | no | -0.0699 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | CDv | 6.642000e-03 | -2.714e-07 | 4.447e-07 | 1.000e-06 | -2.978e-07 | 1.24e-08 | yes (rounding) | -0.0045 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | CMz | -3.954594e-03 | -2.385e-05 | 6.588e-05 | 2.218e-04 | -1.672e-05 | 2.76e-06 | yes | 0.4228 |
| SA 0.01-20 NN Drift | CL | 6.550247e-01 | -5.452e-05 | 1.358e-04 | 4.609e-04 | -5.019e-05 | 5.38e-06 | yes | -0.0077 |
| SA 0.01-20 NN Drift | CD | 8.427870e-03 | 1.057e-05 | 3.067e-05 | 1.170e-04 | -1.655e-06 | 1.34e-06 | no | -0.0196 |
| SA 0.01-20 NN Drift | CDp | 1.785618e-03 | 1.099e-05 | 3.063e-05 | 1.180e-04 | -1.242e-06 | 1.34e-06 | no | -0.0696 |
| SA 0.01-20 NN Drift | CDv | 6.642000e-03 | -2.714e-07 | 4.447e-07 | 1.000e-06 | -2.979e-07 | 1.24e-08 | yes (rounding) | -0.0045 |
| SA 0.01-20 NN Drift | CMz | -3.954594e-03 | -2.483e-05 | 6.310e-05 | 1.918e-04 | -1.688e-05 | 2.63e-06 | yes | 0.4268 |
| SA 0.01-20 NN Drift Seed2 | CL | 6.550247e-01 | -6.475e-05 | 1.494e-04 | 4.175e-04 | -6.598e-05 | 5.63e-06 | yes | -0.0101 |
| SA 0.01-20 NN Drift Seed2 | CD | 8.427870e-03 | 6.319e-06 | 2.968e-05 | 1.219e-04 | -2.197e-06 | 1.29e-06 | no | -0.0261 |
| SA 0.01-20 NN Drift Seed2 | CDp | 1.785618e-03 | 6.819e-06 | 2.967e-05 | 1.220e-04 | -1.762e-06 | 1.29e-06 | no | -0.0987 |
| SA 0.01-20 NN Drift Seed2 | CDv | 6.642000e-03 | -2.814e-07 | 4.497e-07 | 1.000e-06 | -3.047e-07 | 1.23e-08 | yes (rounding) | -0.0046 |
| SA 0.01-20 NN Drift Seed2 | CMz | -3.954594e-03 | -2.476e-05 | 6.672e-05 | 2.213e-04 | -1.848e-05 | 2.77e-06 | yes | 0.4674 |
| SST 0.01-20 ConsGalerkinProj Drift | CL | 6.530779e-01 | -3.829e-06 | 1.544e-04 | 4.713e-04 | -8.140e-06 | 6.74e-06 | no | -0.0012 |
| SST 0.01-20 ConsGalerkinProj Drift | CD | 8.116879e-03 | 1.308e-05 | 3.168e-05 | 1.074e-04 | 1.031e-07 | 1.39e-06 | no | 0.0013 |
| SST 0.01-20 ConsGalerkinProj Drift | CDp | 1.747427e-03 | 1.344e-05 | 3.164e-05 | 1.070e-04 | 3.909e-07 | 1.38e-06 | no | 0.0224 |
| SST 0.01-20 ConsGalerkinProj Drift | CDv | 6.369420e-03 | -4.322e-07 | 5.201e-07 | 2.000e-06 | -3.840e-07 | 1.20e-08 | yes (rounding) | -0.0060 |
| SST 0.01-20 ConsGalerkinProj Drift | CMz | -4.526376e-03 | -1.068e-05 | 8.061e-05 | 2.239e-04 | -5.936e-06 | 3.51e-06 | no | 0.1311 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | CL | 6.530779e-01 | -1.004e-05 | 1.580e-04 | 4.774e-04 | -1.697e-05 | 6.86e-06 | yes | -0.0026 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | CD | 8.116879e-03 | 9.043e-06 | 3.129e-05 | 1.137e-04 | 6.186e-07 | 1.37e-06 | no | 0.0076 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | CDp | 1.747427e-03 | 9.420e-06 | 3.129e-05 | 1.140e-04 | 9.548e-07 | 1.37e-06 | no | 0.0546 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | CDv | 6.369420e-03 | -5.025e-07 | 6.448e-07 | 2.000e-06 | -4.804e-07 | 1.46e-08 | yes (rounding) | -0.0075 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | CMz | -4.526376e-03 | -1.060e-05 | 8.074e-05 | 2.350e-04 | -5.931e-06 | 3.52e-06 | no | 0.1310 |
| SST 0.01-20 NN Drift | CL | 6.530779e-01 | -9.734e-06 | 1.505e-04 | 4.154e-04 | -1.164e-05 | 6.55e-06 | no | -0.0018 |
| SST 0.01-20 NN Drift | CD | 8.116879e-03 | 1.323e-05 | 3.157e-05 | 1.124e-04 | 1.652e-07 | 1.38e-06 | no | 0.0020 |
| SST 0.01-20 NN Drift | CDp | 1.747427e-03 | 1.359e-05 | 3.153e-05 | 1.130e-04 | 4.720e-07 | 1.38e-06 | no | 0.0270 |
| SST 0.01-20 NN Drift | CDv | 6.369420e-03 | -4.422e-07 | 5.403e-07 | 2.000e-06 | -3.985e-07 | 1.25e-08 | yes (rounding) | -0.0063 |
| SST 0.01-20 NN Drift | CMz | -4.526376e-03 | -1.134e-05 | 8.081e-05 | 2.230e-04 | -6.940e-06 | 3.52e-06 | no | 0.1533 |
| SST 0.01-20 NN Drift Seed2 | CL | 6.530779e-01 | -2.057e-05 | 1.548e-04 | 4.576e-04 | -2.568e-05 | 6.65e-06 | yes | -0.0039 |
| SST 0.01-20 NN Drift Seed2 | CD | 8.116879e-03 | 9.219e-06 | 3.112e-05 | 1.207e-04 | 1.805e-07 | 1.36e-06 | no | 0.0022 |
| SST 0.01-20 NN Drift Seed2 | CDp | 1.747427e-03 | 9.623e-06 | 3.113e-05 | 1.210e-04 | 5.448e-07 | 1.36e-06 | no | 0.0312 |
| SST 0.01-20 NN Drift Seed2 | CDv | 6.369420e-03 | -5.050e-07 | 6.487e-07 | 2.000e-06 | -4.835e-07 | 1.47e-08 | yes (rounding) | -0.0076 |
| SST 0.01-20 NN Drift Seed2 | CMz | -4.526376e-03 | -1.191e-05 | 7.876e-05 | 2.321e-04 | -7.965e-06 | 3.42e-06 | yes | 0.1760 |

Drag split at full precision (surface integration of Cp and Cf, `rans_surface_forces.csv`), same OLS trend per 100 steps:

| case | steps | dCL mean | dCL max | dCL trend | dCDp mean | dCDp max | dCDp trend | dCDv mean | dCDv max | dCDv trend |
|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | 478 | -4.708e-05 | 4.482e-04 | -4.730e-05 ± 5.1e-06 * | 1.103e-05 | 1.024e-04 | -1.072e-06 ± 1.2e-06 | -4.414e-07 | 1.376e-06 | -4.113e-07 ± 6.2e-09 * |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | 478 | -5.311e-05 | 4.038e-04 | -5.763e-05 ± 5.2e-06 * | 7.146e-06 | 1.125e-04 | -1.055e-06 ± 1.1e-06 | -4.942e-07 | 1.485e-06 | -4.423e-07 ± 5.9e-09 * |
| SA 0.01-20 NN Drift | 478 | -5.167e-05 | 4.609e-04 | -5.019e-05 ± 4.9e-06 * | 1.156e-05 | 1.183e-04 | -9.488e-07 ± 1.2e-06 | -4.397e-07 | 1.388e-06 | -4.121e-07 ± 6.3e-09 * |
| SA 0.01-20 NN Drift Seed2 | 478 | -6.202e-05 | 4.175e-04 | -6.473e-05 ± 5.2e-06 * | 7.505e-06 | 1.220e-04 | -1.549e-06 ± 1.2e-06 | -5.014e-07 | 1.529e-06 | -4.525e-07 ± 6.1e-09 * |
| SST 0.01-20 ConsGalerkinProj Drift | 478 | -1.840e-07 | 4.713e-04 | -7.650e-06 ± 6.1e-06 | 1.420e-05 | 1.075e-04 | 6.316e-07 ± 1.3e-06 | -3.416e-07 | 1.085e-06 | -2.865e-07 ± 4.2e-09 * |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | 478 | -6.274e-06 | 4.774e-04 | -1.561e-05 ± 6.2e-06 * | 1.044e-05 | 1.138e-04 | 1.379e-06 ± 1.2e-06 | -3.747e-07 | 1.237e-06 | -3.235e-07 ± 4.7e-09 * |
| SST 0.01-20 NN Drift | 478 | -5.477e-06 | 4.154e-04 | -1.141e-05 ± 6.0e-06 | 1.474e-05 | 1.134e-04 | 6.918e-07 ± 1.2e-06 | -3.424e-07 | 1.093e-06 | -2.868e-07 ± 4.2e-09 * |
| SST 0.01-20 NN Drift Seed2 | 478 | -1.668e-05 | 4.576e-04 | -2.372e-05 ± 6.1e-06 * | 1.087e-05 | 1.209e-04 | 1.027e-06 ± 1.2e-06 | -3.785e-07 | 1.263e-06 | -3.293e-07 ± 4.9e-09 * |

(* = resolved at 2 standard errors)

## III.4 Response to a single transfer

One transfer at the event step, then as many steps again. Onset lag = first step after the event at which |deviation| reaches the given fraction of that case's own peak. "pre-floor" is the largest deviation before the event (bit-identity check); "leak" the deviation at the event step. CENSORED = the peak lies in the last 10% of the window.

**CL**

| case | d | event | window | pre-floor | leak | peak | peak lag | onset 1% | onset 5% | onset 10% | onset 25% | onset 50% | predicted 1000·d | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 2.85e-06 | 2.697e-05 | 199 | 0 | 0 | 0 | 2 | 3 | 50 | CENSORED |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 2.53e-07 | 6.854e-05 | 191 | 1 | 1 | 2 | 6 | 9 | 100 | CENSORED |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 1.33e-07 | 2.539e-05 | 191 | 1 | 1 | 2 | 6 | 8 | 100 | CENSORED |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 3.18e-08 | 8.869e-06 | 5 | 1 | 1 | 1 | 2 | 3 | 200 |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 1.26e-08 | 1.423e-05 | 10 | 1 | 1 | 1 | 2 | 5 | 200 |  |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 2.84e-06 | 2.629e-05 | 199 | 0 | 0 | 0 | 2 | 3 | 50 | CENSORED |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 2.52e-07 | 6.746e-05 | 191 | 1 | 1 | 2 | 6 | 9 | 100 | CENSORED |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 1.33e-07 | 2.528e-05 | 191 | 1 | 1 | 2 | 6 | 8 | 100 | CENSORED |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 3.21e-08 | 8.309e-06 | 5 | 1 | 1 | 1 | 2 | 3 | 200 |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 1.28e-08 | 1.384e-05 | 10 | 1 | 1 | 1 | 2 | 5 | 200 |  |

**CDp**

| case | d | event | window | pre-floor | leak | peak | peak lag | onset 1% | onset 5% | onset 10% | onset 25% | onset 50% | predicted 1000·d | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 0.00e+00 | 5.000e-06 | 30 | 1 | 1 | 1 | 1 | 2 | 50 | rounding |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 7.000e-06 | 2 | 1 | 1 | 1 | 1 | 1 | 100 | rounding |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 3.000e-06 | 2 | 1 | 1 | 1 | 1 | 1 | 100 | rounding |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 0.00e+00 | 1.000e-06 | 5 | 2 | 2 | 2 | 2 | 2 | 200 | rounding |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 0.00e+00 | 1.000e-06 | 13 | 2 | 2 | 2 | 2 | 2 | 200 | rounding |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 0.00e+00 | 5.000e-06 | 2 | 1 | 1 | 1 | 1 | 2 | 50 | rounding |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 7.000e-06 | 2 | 1 | 1 | 1 | 1 | 1 | 100 | rounding |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 3.000e-06 | 2 | 1 | 1 | 1 | 1 | 2 | 100 | rounding |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 0.00e+00 | 1.000e-06 | 5 | 2 | 2 | 2 | 2 | 2 | 200 | rounding |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 0.00e+00 | 1.000e-06 | 11 | 2 | 2 | 2 | 2 | 2 | 200 | rounding |

**CDv**

| case | d | event | window | pre-floor | leak | peak | peak lag | onset 1% | onset 5% | onset 10% | onset 25% | onset 50% | predicted 1000·d | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 50 | rounding |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 100 | rounding |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 100 | rounding |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 200 | rounding |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 200 | rounding |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 0.00e+00 | 0.00e+00 | 1.000e-06 | 19 | 19 | 19 | 19 | 19 | 19 | 50 | rounding |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 1.000e-06 | 50 | 50 | 50 | 50 | 50 | 50 | 100 | rounding |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 0.00e+00 | 0.00e+00 | 1.000e-06 | 50 | 50 | 50 | 50 | 50 | 50 | 100 | rounding |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 200 | rounding |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 0.00e+00 | 0.00e+00 | 0.000e+00 | 0 | - | - | - | - | - | 200 | rounding |

**CD**

| case | d | event | window | pre-floor | leak | peak | peak lag | onset 1% | onset 5% | onset 10% | onset 25% | onset 50% | predicted 1000·d | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 1.00e-12 | 2.69e-07 | 4.662e-06 | 5 | 0 | 0 | 1 | 1 | 2 | 50 |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 1.00e-12 | 1.40e-07 | 7.423e-06 | 2 | 0 | 1 | 1 | 1 | 1 | 100 |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 1.00e-12 | 4.93e-08 | 3.258e-06 | 2 | 0 | 1 | 1 | 1 | 1 | 100 |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 1.00e-12 | 5.47e-09 | 9.777e-07 | 31 | 1 | 1 | 1 | 1 | 2 | 200 |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 1.00e-12 | 1.54e-09 | 1.492e-06 | 2 | 1 | 1 | 1 | 2 | 2 | 200 |  |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 1.00e-12 | 2.66e-07 | 4.539e-06 | 5 | 0 | 0 | 1 | 1 | 2 | 50 |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 1.00e-12 | 1.40e-07 | 7.239e-06 | 2 | 0 | 1 | 1 | 1 | 1 | 100 |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 1.00e-12 | 4.92e-08 | 3.196e-06 | 2 | 0 | 1 | 1 | 1 | 1 | 100 |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 1.00e-12 | 5.47e-09 | 9.738e-07 | 31 | 1 | 1 | 1 | 1 | 2 | 200 |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 1.00e-12 | 1.54e-09 | 1.498e-06 | 2 | 1 | 1 | 1 | 2 | 2 | 200 |  |

**CMz**

| case | d | event | window | pre-floor | leak | peak | peak lag | onset 1% | onset 5% | onset 10% | onset 25% | onset 50% | predicted 1000·d | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 1.00e-12 | 5.50e-07 | 7.795e-06 | 42 | 0 | 0 | 1 | 1 | 2 | 50 |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 1.00e-12 | 3.77e-08 | 2.307e-05 | 3 | 1 | 1 | 2 | 2 | 2 | 100 |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 1.00e-12 | 2.24e-08 | 8.490e-06 | 3 | 1 | 2 | 2 | 2 | 2 | 100 |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 1.00e-12 | 3.84e-09 | 2.648e-06 | 4 | 1 | 1 | 3 | 3 | 3 | 200 |  |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 1.00e-12 | 1.22e-09 | 4.891e-06 | 42 | 1 | 2 | 3 | 3 | 4 | 200 |  |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 200 | 199 | 1.00e-12 | 5.55e-07 | 7.425e-06 | 42 | 0 | 0 | 1 | 1 | 2 | 50 |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 200 | 199 | 1.00e-12 | 3.73e-08 | 2.274e-05 | 3 | 1 | 1 | 2 | 2 | 2 | 100 |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 200 | 199 | 1.00e-12 | 2.22e-08 | 8.404e-06 | 3 | 1 | 2 | 2 | 2 | 2 | 100 |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 400 | 399 | 1.00e-12 | 3.91e-09 | 2.382e-06 | 4 | 1 | 1 | 3 | 3 | 3 | 200 |  |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 200 | 199 | 1.00e-12 | 1.29e-09 | 4.759e-06 | 42 | 1 | 2 | 3 | 3 | 4 | 200 |  |

Lag scaling fit lag = C·d + b over the event-at-200 cases (prediction C = 1000):

| model | coeff | threshold | cases | C | b | R² | C / 1000 |
|---|---|---|---|---|---|---|---|
| SA | CL | 0.01 | 3 | 5.7 | 0.00 | 0.5714 | 0.006 |
| SA | CL | 0.05 | 3 | 5.7 | 0.00 | 0.5714 | 0.006 |
| SA | CL | 0.1 | 3 | 4.3 | 0.50 | 0.1071 | 0.004 |
| SA | CL | 0.25 | 3 | -5.7 | 4.00 | 0.0357 | -0.006 |
| SA | CL | 0.5 | 3 | 7.1 | 4.50 | 0.0470 | 0.007 |
| SA | CDp | 0.01 | 3 | 7.1 | 0.50 | 0.8929 | 0.007 |
| SA | CDp | 0.05 | 3 | 7.1 | 0.50 | 0.8929 | 0.007 |
| SA | CDp | 0.1 | 3 | 7.1 | 0.50 | 0.8929 | 0.007 |
| SA | CDp | 0.25 | 3 | 7.1 | 0.50 | 0.8929 | 0.007 |
| SA | CDp | 0.5 | 3 | 1.4 | 1.50 | 0.0357 | 0.001 |
| SST | CL | 0.01 | 3 | 5.7 | 0.00 | 0.5714 | 0.006 |
| SST | CL | 0.05 | 3 | 5.7 | 0.00 | 0.5714 | 0.006 |
| SST | CL | 0.1 | 3 | 4.3 | 0.50 | 0.1071 | 0.004 |
| SST | CL | 0.25 | 3 | -5.7 | 4.00 | 0.0357 | -0.006 |
| SST | CL | 0.5 | 3 | 7.1 | 4.50 | 0.0470 | 0.007 |
| SST | CDp | 0.01 | 3 | 7.1 | 0.50 | 0.8929 | 0.007 |
| SST | CDp | 0.05 | 3 | 7.1 | 0.50 | 0.8929 | 0.007 |
| SST | CDp | 0.1 | 3 | 7.1 | 0.50 | 0.8929 | 0.007 |
| SST | CDp | 0.25 | 3 | 7.1 | 0.50 | 0.8929 | 0.007 |
| SST | CDp | 0.5 | 3 | -0.0 | 2.00 | nan | -0.000 |

**SA amplitude control** (d = 0.10, SafeFactor 2 against 4): CL peak ratio 2.70; onset shift in steps by threshold 1%: +0, 5%: +0, 10%: +0, 25%: +0, 50%: +1.

**SST amplitude control** (d = 0.10, SafeFactor 2 against 4): CL peak ratio 2.67; onset shift in steps by threshold 1%: +0, 5%: +0, 10%: +0, 25%: +0, 50%: +1.

Peak deviations after the event from the surface integration (full-precision split):

| case | dCL pre-floor | dCL peak | at lag | dCDp pre-floor | dCDp peak | at lag | dCDv pre-floor | dCDv peak | at lag |
|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 4.100e-13 | 2.697e-05 | 199 | 1.701e-14 | 4.661e-06 | 5 | 0.000e+00 | 1.679e-08 | 148 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 4.100e-13 | -6.854e-05 | 191 | 1.701e-14 | -7.422e-06 | 2 | 0.000e+00 | 4.767e-08 | 199 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 4.100e-13 | -2.539e-05 | 191 | 1.701e-14 | -3.257e-06 | 2 | 0.000e+00 | 1.796e-08 | 199 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 7.546e-12 | 8.866e-06 | 5 | 6.599e-13 | -9.783e-07 | 31 | 0.000e+00 | -7.508e-09 | 153 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 4.100e-13 | 1.423e-05 | 10 | 1.701e-14 | -1.490e-06 | 2 | 0.000e+00 | -4.284e-09 | 199 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 9.445e-12 | 2.629e-05 | 199 | 3.848e-14 | 4.539e-06 | 5 | 0.000e+00 | 1.685e-08 | 198 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 9.445e-12 | -6.747e-05 | 191 | 3.848e-14 | -7.237e-06 | 2 | 0.000e+00 | 4.649e-08 | 199 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 9.445e-12 | -2.528e-05 | 191 | 3.848e-14 | -3.196e-06 | 2 | 0.000e+00 | 1.589e-08 | 199 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 9.445e-12 | 8.311e-06 | 5 | 1.246e-13 | -9.734e-07 | 31 | 0.000e+00 | -9.273e-09 | 204 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 9.445e-12 | 1.384e-05 | 10 | 3.848e-14 | -1.497e-06 | 2 | 0.000e+00 | -1.241e-09 | 199 |

## III.4b Pressure probe inside the launch annulus

Shell-RMS of the pressure difference case − reference on shells of wall distance r inside the launch radius d (nodes that never moved), as in Part I §6. Onset lag per shell at several fractions of the shell's own peak; the acoustic arrival at a = 1 is 1000·(d − r) steps (`rans_pressure_probe.csv`, `RANS_<MODEL>_PressureProbe_Pressure.png`).

The volume fields are written in single precision; one ULP of the pressure (p ≈ 0.714) is 5.96e-08. Shells whose peak is below 3 ULP are marked: their onset lags are rounding, not arrivals. The force histories are double precision and are not affected.

| case | d | shell r | nodes | d − r | acoustic lag | peak | onset 2% | onset 5% | onset 10% | onset 25% | onset 50% | onset 75% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.005 | 4444 | 0.045 | 45 | 3.392e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.01 | 1556 | 0.040 | 40 | 3.332e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.02 | 839 | 0.030 | 30 | 5.378e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.04 | 548 | 0.010 | 10 | 7.989e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.005 | 4444 | 0.095 | 95 | 3.906e-07 | 0 | 0 | 0 | 1 | 1 | 1 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.01 | 1556 | 0.090 | 90 | 4.902e-07 | 0 | 0 | 0 | 1 | 1 | 1 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.02 | 839 | 0.080 | 80 | 5.238e-07 | 0 | 0 | 0 | 0 | 1 | 1 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.04 | 548 | 0.060 | 60 | 4.023e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.06 | 418 | 0.040 | 40 | 6.139e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.08 | 335 | 0.020 | 20 | 5.678e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.005 | 4444 | 0.095 | 95 | 1.307e-07 (< 3 ULP) | 0 | 0 | 0 | 1 | 1 | 1 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.01 | 1556 | 0.090 | 90 | 1.682e-07 (< 3 ULP) | 0 | 0 | 0 | 1 | 1 | 1 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.02 | 839 | 0.080 | 80 | 1.762e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 1 | 1 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.04 | 548 | 0.060 | 60 | 1.403e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.06 | 418 | 0.040 | 40 | 2.102e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.08 | 335 | 0.020 | 20 | 2.321e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.005 | 4444 | 0.195 | 195 | 7.433e-08 (< 3 ULP) | 0 | 0 | 0 | 1 | 2 | 26 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.01 | 1556 | 0.190 | 190 | 6.136e-08 (< 3 ULP) | 0 | 0 | 0 | 1 | 1 | 23 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.02 | 839 | 0.180 | 180 | 4.956e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 1 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.04 | 548 | 0.160 | 160 | 4.253e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.06 | 418 | 0.140 | 140 | 4.674e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.08 | 335 | 0.120 | 120 | 4.753e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.1 | 228 | 0.100 | 100 | 5.147e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.13 | 162 | 0.070 | 70 | 6.352e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.16 | 109 | 0.040 | 40 | 1.106e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.18 | 121 | 0.020 | 20 | 7.194e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.005 | 4444 | 0.195 | 195 | 6.036e-08 (< 3 ULP) | 0 | 0 | 0 | 1 | 1 | 2 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.01 | 1556 | 0.190 | 190 | 6.040e-08 (< 3 ULP) | 0 | 0 | 0 | 1 | 1 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.02 | 839 | 0.180 | 180 | 6.208e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 1 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.04 | 548 | 0.160 | 160 | 6.355e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 1 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.06 | 418 | 0.140 | 140 | 6.327e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 1 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.08 | 335 | 0.120 | 120 | 6.497e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.1 | 228 | 0.100 | 100 | 8.756e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.13 | 162 | 0.070 | 70 | 1.038e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.16 | 109 | 0.040 | 40 | 1.520e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.18 | 121 | 0.020 | 20 | 6.146e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.005 | 4444 | 0.045 | 45 | 3.389e-07 | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.01 | 1556 | 0.040 | 40 | 3.316e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.02 | 839 | 0.030 | 30 | 5.423e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 0.05 | 0.04 | 548 | 0.010 | 10 | 8.122e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.005 | 4444 | 0.095 | 95 | 3.872e-07 | 0 | 0 | 0 | 1 | 1 | 1 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.01 | 1556 | 0.090 | 90 | 4.768e-07 | 0 | 0 | 0 | 1 | 1 | 1 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.02 | 839 | 0.080 | 80 | 5.096e-07 | 0 | 0 | 0 | 0 | 1 | 1 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.04 | 548 | 0.060 | 60 | 3.922e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.06 | 418 | 0.040 | 40 | 5.975e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 0.1 | 0.08 | 335 | 0.020 | 20 | 5.555e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.005 | 4444 | 0.095 | 95 | 1.298e-07 (< 3 ULP) | 0 | 0 | 0 | 1 | 1 | 1 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.01 | 1556 | 0.090 | 90 | 1.651e-07 (< 3 ULP) | 0 | 0 | 0 | 1 | 1 | 1 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.02 | 839 | 0.080 | 80 | 1.731e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 1 | 1 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.04 | 548 | 0.060 | 60 | 1.379e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.06 | 418 | 0.040 | 40 | 2.067e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 0.1 | 0.08 | 335 | 0.020 | 20 | 2.340e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.005 | 4444 | 0.195 | 195 | 7.387e-08 (< 3 ULP) | 0 | 0 | 0 | 1 | 2 | 26 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.01 | 1556 | 0.190 | 190 | 5.905e-08 (< 3 ULP) | 0 | 0 | 0 | 1 | 1 | 23 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.02 | 839 | 0.180 | 180 | 5.045e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 1 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.04 | 548 | 0.160 | 160 | 4.291e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.06 | 418 | 0.140 | 140 | 4.235e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.08 | 335 | 0.120 | 120 | 4.685e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.1 | 228 | 0.100 | 100 | 4.818e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.13 | 162 | 0.070 | 70 | 6.034e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.16 | 109 | 0.040 | 40 | 1.048e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 0.2 | 0.18 | 121 | 0.020 | 20 | 7.010e-07 | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.005 | 4444 | 0.195 | 195 | 6.015e-08 (< 3 ULP) | 0 | 0 | 0 | 1 | 1 | 2 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.01 | 1556 | 0.190 | 190 | 5.993e-08 (< 3 ULP) | 0 | 0 | 0 | 1 | 1 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.02 | 839 | 0.180 | 180 | 6.125e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 1 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.04 | 548 | 0.160 | 160 | 5.977e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 1 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.06 | 418 | 0.140 | 140 | 6.473e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 1 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.08 | 335 | 0.120 | 120 | 6.626e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.1 | 228 | 0.100 | 100 | 8.233e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.13 | 162 | 0.070 | 70 | 9.540e-08 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.16 | 109 | 0.040 | 40 | 1.425e-07 (< 3 ULP) | 0 | 0 | 0 | 0 | 0 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 0.2 | 0.18 | 121 | 0.020 | 20 | 5.975e-07 | 0 | 0 | 0 | 0 | 0 | 0 |

Differential fit lag = S·(d − r) + b over shells with d − r ≥ 0.08 (prediction S = 1000, speed 1):

| case | threshold | shells | S | b | R² | implied speed |
|---|---|---|---|---|---|---|
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 2% | 3 | all shells at lag 0 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 5% | 3 | all shells at lag 0 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 10% | 3 | all shells at lag 0 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 25% | 3 | 71 | -5.6 | 0.893 | 14.000 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 50% | 3 | all shells at lag 1 |  |  |  |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 75% | 3 | all shells at lag 1 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 2% | 3 | all shells at lag 0 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 5% | 3 | all shells at lag 0 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 10% | 3 | all shells at lag 0 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 25% | 3 | 71 | -5.6 | 0.893 | 14.000 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 50% | 3 | all shells at lag 1 |  |  |  |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 75% | 3 | all shells at lag 1 |  |  |  |

## III.5 Surface pressure and skin friction

RMS and maximum over the 370 airfoil nodes of the difference to the reference at the same step (`rans_surface_norms.csv`). Cf_x is the streamwise wall-shear coefficient. The minimum-Cf_x node of the reference is listed as a separation check (a negative Cf_x at the leading edge is the stagnation region below the nose at AoA 6°, not separation).

| case | when | step | RMS dCf_x | max dCf_x | x of max | RMS dCp | max dCp | x of max | max dy+ | min Cf_x ref | x | max y+ ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | step 9 | 9 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9972e-03 | 0.0018 | 2.694 |
| SA 0.01-20 ConsGalerkinProj Drift | step 200 | 200 | 2.698e-06 | 1.719e-05 | 0.3821 | 2.838e-03 | 2.753e-02 | 0.5220 | 3.828e-03 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.01-20 ConsGalerkinProj Drift | last | 399 | 3.480e-06 | 1.978e-05 | 0.4507 | 2.742e-03 | 1.832e-02 | 0.4580 | 3.649e-03 | -2.9965e-03 | 0.0018 | 2.694 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | step 9 | 9 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9972e-03 | 0.0018 | 2.694 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | step 200 | 200 | 2.998e-06 | 1.807e-05 | 0.4283 | 2.077e-03 | 1.230e-02 | 0.6014 | 3.289e-03 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | last | 399 | 3.132e-06 | 1.500e-05 | 0.4580 | 3.889e-03 | 3.452e-02 | 0.4507 | 2.742e-03 | -2.9965e-03 | 0.0018 | 2.694 |
| SA 0.01-20 NN Drift | step 9 | 9 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9972e-03 | 0.0018 | 2.694 |
| SA 0.01-20 NN Drift | step 200 | 200 | 2.729e-06 | 1.765e-05 | 0.3821 | 2.876e-03 | 2.736e-02 | 0.5220 | 3.869e-03 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.01-20 NN Drift | last | 399 | 3.494e-06 | 1.960e-05 | 0.4507 | 2.777e-03 | 1.815e-02 | 0.4580 | 3.616e-03 | -2.9965e-03 | 0.0018 | 2.694 |
| SA 0.01-20 NN Drift Seed2 | step 9 | 9 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9972e-03 | 0.0018 | 2.694 |
| SA 0.01-20 NN Drift Seed2 | step 200 | 200 | 3.025e-06 | 1.831e-05 | 0.4283 | 2.089e-03 | 1.170e-02 | 0.5151 | 3.303e-03 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.01-20 NN Drift Seed2 | last | 399 | 3.166e-06 | 1.515e-05 | 0.4580 | 3.897e-03 | 3.458e-02 | 0.4507 | 2.770e-03 | -2.9965e-03 | 0.0018 | 2.694 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | event+1 | 201 | 2.263e-09 | 1.444e-08 | 0.0730 | 7.538e-05 | 4.785e-04 | 1.0000 | 6.914e-06 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | event+50 | 250 | 4.518e-08 | 2.682e-07 | 0.0085 | 4.776e-05 | 2.470e-04 | 0.0001 | 8.762e-05 | -2.9968e-03 | 0.0018 | 2.694 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | last | 399 | 1.447e-07 | 1.891e-06 | 0.0085 | 4.690e-05 | 2.022e-04 | 0.0001 | 6.181e-04 | -2.9965e-03 | 0.0018 | 2.694 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | event-1 | 199 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | event+1 | 201 | 9.839e-10 | 5.588e-09 | 0.0564 | 2.601e-05 | 1.363e-04 | 1.0000 | 1.729e-06 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | event+50 | 250 | 4.248e-08 | 6.042e-07 | 0.0085 | 7.255e-05 | 3.932e-04 | 0.0007 | 1.972e-04 | -2.9968e-03 | 0.0018 | 2.694 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | last | 399 | 3.134e-07 | 5.164e-06 | 0.0085 | 9.488e-05 | 4.738e-04 | 0.0007 | 1.689e-03 | -2.9965e-03 | 0.0018 | 2.694 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | event+1 | 201 | 3.945e-10 | 2.328e-09 | 0.0730 | 9.685e-06 | 4.312e-05 | 1.0000 | 5.364e-07 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | event+50 | 250 | 1.522e-08 | 2.200e-07 | 0.0085 | 2.236e-05 | 1.197e-04 | 0.0018 | 7.182e-05 | -2.9968e-03 | 0.0018 | 2.694 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | last | 399 | 1.132e-07 | 1.934e-06 | 0.0085 | 3.037e-05 | 1.471e-04 | 0.0007 | 6.321e-04 | -2.9965e-03 | 0.0018 | 2.694 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | event-1 | 399 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9965e-03 | 0.0018 | 2.694 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | event+1 | 401 | 1.182e-10 | 4.657e-10 | 0.0012 | 1.943e-06 | 7.153e-06 | 0.1115 | 2.384e-07 | -2.9964e-03 | 0.0018 | 2.694 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | event+50 | 450 | 6.664e-09 | 4.983e-08 | 0.0085 | 1.666e-05 | 8.595e-05 | 0.0001 | 1.627e-05 | -2.9964e-03 | 0.0018 | 2.694 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | last | 799 | 4.949e-08 | 7.919e-07 | 0.0085 | 1.214e-05 | 6.056e-05 | 0.0001 | 2.589e-04 | -2.9961e-03 | 0.0018 | 2.693 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | event+1 | 201 | 1.831e-10 | 9.313e-10 | 0.0051 | 2.210e-06 | 9.298e-06 | 0.0429 | 2.384e-07 | -2.9969e-03 | 0.0018 | 2.694 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | event+50 | 250 | 3.569e-09 | 2.421e-08 | 1.0000 | 7.308e-06 | 2.491e-05 | 0.0007 | 6.974e-06 | -2.9968e-03 | 0.0018 | 2.694 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | last | 399 | 1.009e-08 | 5.320e-08 | 0.0085 | 9.190e-06 | 3.123e-05 | 0.0273 | 1.740e-05 | -2.9965e-03 | 0.0018 | 2.694 |
| SST 0.01-20 ConsGalerkinProj Drift | step 9 | 9 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9245e-03 | 0.0018 | 2.642 |
| SST 0.01-20 ConsGalerkinProj Drift | step 200 | 200 | 2.647e-06 | 1.408e-05 | 0.5220 | 3.149e-03 | 3.058e-02 | 0.5220 | 3.985e-03 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.01-20 ConsGalerkinProj Drift | last | 399 | 2.806e-06 | 1.770e-05 | 0.4507 | 2.768e-03 | 1.806e-02 | 0.4580 | 3.329e-03 | -2.9240e-03 | 0.0018 | 2.642 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | step 9 | 9 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9245e-03 | 0.0018 | 2.642 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | step 200 | 200 | 2.727e-06 | 1.515e-05 | 0.4283 | 2.174e-03 | 1.633e-02 | 0.6014 | 2.987e-03 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | last | 399 | 2.774e-06 | 1.246e-05 | 0.4580 | 3.360e-03 | 2.453e-02 | 0.4580 | 2.465e-03 | -2.9240e-03 | 0.0018 | 2.642 |
| SST 0.01-20 NN Drift | step 9 | 9 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9245e-03 | 0.0018 | 2.642 |
| SST 0.01-20 NN Drift | step 200 | 200 | 2.671e-06 | 1.408e-05 | 0.5220 | 3.182e-03 | 3.029e-02 | 0.5220 | 4.030e-03 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.01-20 NN Drift | last | 399 | 2.823e-06 | 1.776e-05 | 0.4507 | 2.789e-03 | 1.784e-02 | 0.4580 | 3.341e-03 | -2.9240e-03 | 0.0018 | 2.642 |
| SST 0.01-20 NN Drift Seed2 | step 9 | 9 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9245e-03 | 0.0018 | 2.642 |
| SST 0.01-20 NN Drift Seed2 | step 200 | 200 | 2.765e-06 | 1.530e-05 | 0.4283 | 2.158e-03 | 1.589e-02 | 0.6014 | 3.016e-03 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.01-20 NN Drift Seed2 | last | 399 | 2.815e-06 | 1.268e-05 | 0.4580 | 3.379e-03 | 2.459e-02 | 0.4580 | 2.481e-03 | -2.9240e-03 | 0.0018 | 2.642 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | event+1 | 201 | 2.225e-09 | 1.397e-08 | 0.0730 | 7.625e-05 | 4.838e-04 | 1.0000 | 8.225e-06 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | event+50 | 250 | 3.997e-08 | 2.758e-07 | 0.0085 | 4.781e-05 | 2.478e-04 | 0.0001 | 9.161e-05 | -2.9242e-03 | 0.0018 | 2.642 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | last | 399 | 1.618e-07 | 1.940e-06 | 0.0085 | 4.696e-05 | 2.046e-04 | 0.0001 | 6.451e-04 | -2.9240e-03 | 0.0018 | 2.642 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | event-1 | 199 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | event+1 | 201 | 9.620e-10 | 5.588e-09 | 0.0730 | 2.572e-05 | 1.321e-04 | 1.0000 | 1.907e-06 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | event+50 | 250 | 4.186e-08 | 6.200e-07 | 0.0085 | 7.220e-05 | 3.917e-04 | 0.0007 | 2.059e-04 | -2.9242e-03 | 0.0018 | 2.642 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | last | 399 | 3.242e-07 | 5.269e-06 | 0.0085 | 9.427e-05 | 4.727e-04 | 0.0007 | 1.753e-03 | -2.9240e-03 | 0.0018 | 2.642 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | event+1 | 201 | 4.109e-10 | 2.328e-09 | 0.0730 | 9.633e-06 | 4.256e-05 | 0.0730 | 5.960e-07 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | event+50 | 250 | 1.482e-08 | 2.258e-07 | 0.0085 | 2.233e-05 | 1.200e-04 | 0.0018 | 7.498e-05 | -2.9242e-03 | 0.0018 | 2.642 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | last | 399 | 1.151e-07 | 1.975e-06 | 0.0085 | 3.028e-05 | 1.470e-04 | 0.0007 | 6.566e-04 | -2.9240e-03 | 0.0018 | 2.642 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | event-1 | 399 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9240e-03 | 0.0018 | 2.642 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | event+1 | 401 | 1.326e-10 | 9.313e-10 | 0.0067 | 1.938e-06 | 7.167e-06 | 0.1115 | 1.192e-07 | -2.9240e-03 | 0.0018 | 2.642 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | event+50 | 450 | 5.996e-09 | 5.076e-08 | 0.0085 | 1.648e-05 | 8.500e-05 | 0.0001 | 1.687e-05 | -2.9239e-03 | 0.0018 | 2.642 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | last | 799 | 4.864e-08 | 8.029e-07 | 0.0085 | 1.198e-05 | 5.829e-05 | 0.0000 | 2.670e-04 | -2.9236e-03 | 0.0018 | 2.642 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | event-1 | 199 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | 0.000e+00 | 0.0000 | 0.000e+00 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | event+1 | 201 | 1.841e-10 | 9.313e-10 | 0.0037 | 2.212e-06 | 9.298e-06 | 0.0429 | 2.384e-07 | -2.9243e-03 | 0.0018 | 2.642 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | event+50 | 250 | 3.114e-09 | 2.421e-08 | 1.0000 | 7.283e-06 | 2.527e-05 | 0.0007 | 6.497e-06 | -2.9242e-03 | 0.0018 | 2.642 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | last | 399 | 1.079e-08 | 6.147e-08 | 0.0194 | 8.965e-06 | 3.076e-05 | 0.0273 | 1.860e-05 | -2.9240e-03 | 0.0018 | 2.642 |

Distributions (Cf_x and Cp of case and reference at the key steps) are in `rans_surface_profiles.csv` and plotted in `Figures/RANS/<case>_Cf_Cp.png`.

## III.6 Turbulence state

**Transfer defects of the transported variables.** For every transfer, the domain integral of the rho-weighted variable on the source mesh and on the target mesh (`rans_transfer_integrals.csv`, from the interpolation's own `Integrals_<step>.csv`). Only the newest state of each event is summarised (the BDF2 level is transferred too and behaves the same). The last column is the change of the source integral between consecutive events: it is what the transport equations did in the 10 steps in between, and it is listed so that it is not mistaken for a transfer defect.

| case | field | events | max |defect| | median |defect| | min change between events | max change between events | net change first→last event |
|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | Density | 40 | 2.2e-13 | 7.2e-14 | -3.823e-11 | 4.886e-11 | 1.342e-10 |
| SA 0.01-20 ConsGalerkinProj Drift | Energy | 40 | 2.9e-13 | 1.2e-13 | -6.197e-11 | 6.016e-11 | -8.208e-11 |
| SA 0.01-20 ConsGalerkinProj Drift | Nu_Tilde | 40 | 4.1e-09 | 1.2e-09 | -7.731e-07 | 1.985e-07 | -1.660e-05 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | Density | 40 | 1.7e-13 | 5.5e-14 | -3.846e-11 | 7.893e-11 | 2.968e-10 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | Energy | 40 | 2.5e-13 | 7.8e-14 | -6.074e-11 | 1.007e-10 | 1.467e-10 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | Nu_Tilde | 40 | 3.9e-09 | 8.3e-10 | -7.271e-07 | 2.285e-07 | -1.480e-05 |
| SA 0.01-20 NN Drift | Density | 40 | 5.4e-11 | 1.3e-11 | -8.413e-11 | 1.413e-10 | 6.606e-10 |
| SA 0.01-20 NN Drift | Energy | 40 | 3.9e-11 | 5.4e-12 | -1.162e-10 | 1.813e-10 | 5.785e-10 |
| SA 0.01-20 NN Drift | Nu_Tilde | 40 | 2.5e-06 | 6.5e-07 | -3.146e-06 | 9.636e-07 | -2.191e-05 |
| SA 0.01-20 NN Drift Seed2 | Density | 40 | 5.5e-11 | 1.7e-11 | -1.360e-10 | 7.526e-11 | 4.751e-11 |
| SA 0.01-20 NN Drift Seed2 | Energy | 40 | 4.3e-11 | 1.4e-11 | -1.555e-10 | 6.786e-11 | -2.017e-10 |
| SA 0.01-20 NN Drift Seed2 | Nu_Tilde | 40 | 1.5e-06 | 5.3e-07 | -2.008e-06 | 1.397e-06 | -1.216e-05 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | Density | 2 | 1.8e-14 | 1.7e-14 | 1.222e-10 | 1.222e-10 | 1.221e-10 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | Energy | 2 | 1.4e-14 | 1.1e-14 | 9.555e-12 | 9.555e-12 | 9.615e-12 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | Nu_Tilde | 2 | 1.4e-13 | 1.1e-13 | 7.121e-07 | 7.121e-07 | 7.121e-07 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | Density | 2 | 3.4e-15 | 2.3e-15 | 1.190e-10 | 1.190e-10 | 1.189e-10 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | Energy | 2 | 3.6e-15 | 1.8e-15 | 7.778e-12 | 7.778e-12 | 7.834e-12 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | Nu_Tilde | 2 | 1.6e-15 | 1.6e-15 | 9.487e-08 | 9.487e-08 | 9.487e-08 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | Density | 2 | 3.4e-15 | 1.9e-15 | 1.201e-10 | 1.201e-10 | 1.202e-10 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | Energy | 2 | 3.8e-15 | 2.8e-15 | 8.862e-12 | 8.862e-12 | 8.903e-12 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | Nu_Tilde | 2 | 5.3e-15 | 4.0e-15 | -1.481e-07 | -1.481e-07 | -1.481e-07 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | Density | 2 | 2.1e-15 | 1.0e-15 | 3.761e-10 | 3.761e-10 | 3.759e-10 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | Energy | 2 | 1.7e-15 | 1.0e-15 | 8.921e-12 | 8.921e-12 | 8.902e-12 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | Nu_Tilde | 2 | 2.6e-15 | 1.4e-15 | 2.302e-07 | 2.302e-07 | 2.302e-07 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | Density | 2 | 2.3e-15 | 1.6e-15 | 1.211e-10 | 1.211e-10 | 1.212e-10 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | Energy | 2 | 2.7e-15 | 1.9e-15 | 1.070e-11 | 1.070e-11 | 1.086e-11 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | Nu_Tilde | 2 | 6.1e-15 | 4.9e-15 | 1.075e-07 | 1.075e-07 | 1.075e-07 |
| SST 0.01-20 ConsGalerkinProj Drift | Density | 40 | 2.3e-13 | 7.4e-14 | -3.838e-11 | 4.252e-11 | -6.951e-11 |
| SST 0.01-20 ConsGalerkinProj Drift | Energy | 40 | 3.0e-13 | 1.1e-13 | -5.373e-11 | 5.855e-11 | -1.111e-10 |
| SST 0.01-20 ConsGalerkinProj Drift | Turb_Kin_Energy | 40 | 2.6e-08 | 1.1e-08 | 2.732e-07 | 1.684e-06 | 3.367e-05 |
| SST 0.01-20 ConsGalerkinProj Drift | Omega | 40 | 1.1e-09 | 3.6e-10 | 1.885e-08 | 6.912e-08 | 1.692e-06 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Density | 40 | 1.7e-13 | 5.5e-14 | -3.979e-11 | 7.456e-11 | 8.577e-11 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Energy | 40 | 2.3e-13 | 7.3e-14 | -5.729e-11 | 1.043e-10 | 1.113e-10 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Turb_Kin_Energy | 40 | 2.5e-08 | 9.4e-09 | 5.095e-08 | 1.552e-06 | 3.563e-05 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Omega | 40 | 9.1e-10 | 3.4e-10 | 1.271e-08 | 6.754e-08 | 1.683e-06 |
| SST 0.01-20 NN Drift | Density | 40 | 4.9e-11 | 1.2e-11 | -8.615e-11 | 1.319e-10 | 4.336e-10 |
| SST 0.01-20 NN Drift | Energy | 40 | 3.8e-11 | 5.4e-12 | -1.173e-10 | 1.792e-10 | 5.265e-10 |
| SST 0.01-20 NN Drift | Turb_Kin_Energy | 40 | 3.0e-06 | 8.0e-07 | -2.511e-06 | 2.416e-06 | 2.511e-05 |
| SST 0.01-20 NN Drift | Omega | 40 | 8.1e-08 | 2.6e-08 | -4.773e-08 | 8.819e-08 | 1.483e-06 |
| SST 0.01-20 NN Drift Seed2 | Density | 40 | 4.8e-11 | 1.6e-11 | -1.307e-10 | 6.819e-11 | -1.562e-10 |
| SST 0.01-20 NN Drift Seed2 | Energy | 40 | 4.2e-11 | 1.4e-11 | -1.452e-10 | 6.879e-11 | -2.279e-10 |
| SST 0.01-20 NN Drift Seed2 | Turb_Kin_Energy | 40 | 1.8e-06 | 7.0e-07 | -7.170e-07 | 2.934e-06 | 3.453e-05 |
| SST 0.01-20 NN Drift Seed2 | Omega | 40 | 5.3e-08 | 1.8e-08 | -1.897e-08 | 1.077e-07 | 1.679e-06 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Density | 2 | 2.0e-14 | 1.8e-14 | 1.689e-12 | 1.689e-12 | 1.594e-12 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Energy | 2 | 1.6e-14 | 1.5e-14 | 9.596e-12 | 9.596e-12 | 9.615e-12 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Turb_Kin_Energy | 2 | 7.8e-14 | 7.0e-14 | 3.985e-06 | 3.985e-06 | 3.985e-06 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Omega | 2 | 4.5e-13 | 3.6e-13 | 4.577e-08 | 4.577e-08 | 4.577e-08 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Density | 2 | 2.3e-15 | 2.0e-15 | -1.492e-12 | -1.492e-12 | -1.594e-12 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Energy | 2 | 3.8e-15 | 2.4e-15 | 7.868e-12 | 7.868e-12 | 7.834e-12 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Turb_Kin_Energy | 2 | 1.3e-10 | 6.5e-11 | 1.518e-06 | 1.518e-06 | 1.518e-06 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Omega | 2 | 1.4e-13 | 1.4e-13 | -1.318e-09 | -1.318e-09 | -1.318e-09 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Density | 2 | 1.2e-15 | 8.1e-16 | -3.545e-13 | -3.545e-13 | -3.189e-13 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Energy | 2 | 4.7e-15 | 2.5e-15 | 8.936e-12 | 8.936e-12 | 8.903e-12 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Turb_Kin_Energy | 2 | 5.7e-14 | 5.6e-14 | 3.755e-07 | 3.755e-07 | 3.755e-07 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Omega | 2 | 1.5e-13 | 1.4e-13 | -1.775e-08 | -1.775e-08 | -1.775e-08 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Density | 2 | 2.1e-15 | 1.4e-15 | -1.279e-11 | -1.279e-11 | -1.307e-11 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Energy | 2 | 2.3e-15 | 1.7e-15 | 8.088e-12 | 8.088e-12 | 8.191e-12 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Turb_Kin_Energy | 2 | 5.7e-14 | 5.6e-14 | 1.021e-06 | 1.021e-06 | 1.021e-06 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Omega | 2 | 1.4e-13 | 1.4e-13 | -1.523e-08 | -1.523e-08 | -1.523e-08 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Density | 2 | 2.3e-15 | 2.0e-15 | 6.528e-13 | 6.528e-13 | 6.377e-13 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Energy | 2 | 1.6e-15 | 1.2e-15 | 1.076e-11 | 1.076e-11 | 1.068e-11 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Turb_Kin_Energy | 2 | 5.9e-14 | 5.7e-14 | 8.790e-07 | 8.790e-07 | 8.790e-07 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Omega | 2 | 1.5e-13 | 1.4e-13 | -5.107e-09 | -5.107e-09 | -5.107e-09 |

**Minimum of the transferred turbulence variables** over all transfers of each case (source mesh and target mesh, from the interpolation's `Minimums_<step>.csv`). These are the rho-weighted values the interpolation checks (rho*nu_tilde, rho*k, rho*omega), so the free-stream floor of k appears as 0.98e-10 rather than 1e-10; the clipping of III.6 acts on the primitive variables after the division by rho, before this check. "files with target < 0" counts the transfers that produced a negative value somewhere before clipping was introduced (SA), or after it (SST, expected 0):

| case | field | transfers | min source | min target | files with target < 0 |
|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | Nu_Tilde | 80 | 0.000e+00 | -1.176e-16 | 70 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | Nu_Tilde | 80 | 0.000e+00 | -7.337e-17 | 71 |
| SA 0.01-20 NN Drift | Nu_Tilde | 80 | 0.000e+00 | 0.000e+00 | 0 |
| SA 0.01-20 NN Drift Seed2 | Nu_Tilde | 80 | 0.000e+00 | 0.000e+00 | 0 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | Nu_Tilde | 4 | 0.000e+00 | -5.900e-18 | 4 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | Nu_Tilde | 4 | 0.000e+00 | -6.435e-18 | 4 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | Nu_Tilde | 4 | 0.000e+00 | -5.948e-18 | 4 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | Nu_Tilde | 4 | 0.000e+00 | -5.837e-18 | 4 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | Nu_Tilde | 4 | 0.000e+00 | -5.825e-18 | 4 |
| SST 0.01-20 ConsGalerkinProj Drift | Turb_Kin_Energy | 80 | 9.818e-11 | 9.818e-11 | 0 |
| SST 0.01-20 ConsGalerkinProj Drift | Omega | 80 | 7.676e-03 | 7.676e-03 | 0 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Turb_Kin_Energy | 80 | 9.818e-11 | 9.818e-11 | 0 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | Omega | 80 | 7.676e-03 | 7.676e-03 | 0 |
| SST 0.01-20 NN Drift | Turb_Kin_Energy | 80 | 9.818e-11 | 9.818e-11 | 0 |
| SST 0.01-20 NN Drift | Omega | 80 | 7.676e-03 | 7.676e-03 | 0 |
| SST 0.01-20 NN Drift Seed2 | Turb_Kin_Energy | 80 | 9.818e-11 | 9.818e-11 | 0 |
| SST 0.01-20 NN Drift Seed2 | Omega | 80 | 7.676e-03 | 7.676e-03 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Turb_Kin_Energy | 4 | 9.818e-11 | 9.818e-11 | 0 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | Omega | 4 | 7.676e-03 | 7.676e-03 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Turb_Kin_Energy | 4 | 9.818e-11 | 9.818e-11 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | Omega | 4 | 7.676e-03 | 7.676e-03 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Turb_Kin_Energy | 4 | 9.818e-11 | 9.818e-11 | 0 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | Omega | 4 | 7.676e-03 | 7.676e-03 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Turb_Kin_Energy | 4 | 9.818e-11 | 9.818e-11 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | Omega | 4 | 7.676e-03 | 7.676e-03 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Turb_Kin_Energy | 4 | 9.818e-11 | 9.818e-11 | 0 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | Omega | 4 | 7.676e-03 | 7.676e-03 | 0 |

**Clipping to the solver's lower limits inside the interpolation** (`clipFieldsToMinimum` / `FieldsMinimum` in InterpSolution_Params.py): after each transfer, values of the listed fields below their floor are raised to it. The floors are the ones SU2 applies to its SST variables (k ≥ 1e-10, ω ≥ 1e-4). Per case and field: number of transfers, transfers with at least one raised node, largest number of raised nodes in one transfer, lowest value seen before clipping, largest total amount added in one transfer. The SA cases ran with the option off (their lowest transferred nu_tilde was of order -1e-16, see the minima table above).

| case | setting | field | transfers | floor | with raises | max nodes | lowest before | max added |
|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | option absent (old code) | - | 80 | - | - | - | - | - |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | option absent (old code) | - | 80 | - | - | - | - | - |
| SA 0.01-20 NN Drift | option absent (old code) | - | 80 | - | - | - | - | - |
| SA 0.01-20 NN Drift Seed2 | option absent (old code) | - | 80 | - | - | - | - | - |
| SA 0.05-0.10 ConsGalerkinProj OneShot | option absent (old code) | - | 4 | - | - | - | - | - |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | option absent (old code) | - | 4 | - | - | - | - | - |
| SA 0.10-0.15 ConsGalerkinProj OneShot | option absent (old code) | - | 4 | - | - | - | - | - |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | option absent (old code) | - | 4 | - | - | - | - | - |
| SA 0.20-0.25 ConsGalerkinProj OneShot | option absent (old code) | - | 4 | - | - | - | - | - |
| SST 0.01-20 ConsGalerkinProj Drift | ['Turb_Kin_Energy', 'Omega'] | Omega | 80 | 0.0001 | 0 | 0 | 7.68e-03 | 0.00e+00 |
| SST 0.01-20 ConsGalerkinProj Drift | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 80 | 1e-10 | 80 | 423 | -1.29e-07 | 5.67e-07 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | ['Turb_Kin_Energy', 'Omega'] | Omega | 80 | 0.0001 | 0 | 0 | 7.68e-03 | 0.00e+00 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 80 | 1e-10 | 80 | 415 | -1.62e-07 | 6.65e-07 |
| SST 0.01-20 NN Drift | ['Turb_Kin_Energy', 'Omega'] | Omega | 80 | 0.0001 | 0 | 0 | 7.68e-03 | 0.00e+00 |
| SST 0.01-20 NN Drift | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 80 | 1e-10 | 80 | 19 | -2.38e-07 | 5.22e-07 |
| SST 0.01-20 NN Drift Seed2 | ['Turb_Kin_Energy', 'Omega'] | Omega | 80 | 0.0001 | 0 | 0 | 7.68e-03 | 0.00e+00 |
| SST 0.01-20 NN Drift Seed2 | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 80 | 1e-10 | 80 | 13 | -2.17e-07 | 3.82e-07 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Omega | 4 | 0.0001 | 0 | 0 | 7.68e-03 | 0.00e+00 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 4 | 1e-10 | 4 | 213 | 1.00e-10 | 1.17e-15 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | ['Turb_Kin_Energy', 'Omega'] | Omega | 4 | 0.0001 | 0 | 0 | 7.68e-03 | 0.00e+00 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 4 | 1e-10 | 4 | 214 | -2.55e-09 | 2.65e-09 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Omega | 4 | 0.0001 | 0 | 0 | 7.68e-03 | 0.00e+00 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 4 | 1e-10 | 4 | 213 | 1.00e-10 | 1.17e-15 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | ['Turb_Kin_Energy', 'Omega'] | Omega | 4 | 0.0001 | 0 | 0 | 7.68e-03 | 0.00e+00 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 4 | 1e-10 | 4 | 213 | 1.00e-10 | 1.17e-15 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Omega | 4 | 0.0001 | 0 | 0 | 7.68e-03 | 0.00e+00 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | ['Turb_Kin_Energy', 'Omega'] | Turb_Kin_Energy | 4 | 1e-10 | 4 | 213 | 1.00e-10 | 1.17e-15 |

**SA: node-wise discrepancy on never-displaced nodes** (`rans_turb_state_norms_SA.csv`). Nodes whose coordinates coincide with the reference to 1e-12 are compared directly; "band" = wall distance < 0.01 (the wall-layer nodes), "outer" = the other fixed nodes. Values are RMS and max of |case − reference| divided by the reference maximum over the same node set.

| case | step | fixed nodes | of which band | moved nodes | Eddy_Viscosity band RMS | band max | outer RMS | outer max | Nu_Tilde band RMS | band max | outer RMS | outer max |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SA 0.01-20 ConsGalerkinProj Drift | 10 | 14275 | 12233 | 15681 | 1.295e-04 | 3.465e-03 | 7.510e-05 | 3.110e-03 | 1.340e-04 | 3.768e-03 | 2.656e-05 | 9.655e-04 |
| SA 0.01-20 ConsGalerkinProj Drift | 200 | 14275 | 12233 | 15681 | 4.458e-03 | 3.248e-02 | 3.486e-04 | 9.482e-03 | 4.457e-03 | 3.246e-02 | 1.645e-04 | 5.139e-03 |
| SA 0.01-20 ConsGalerkinProj Drift | 395 | 14275 | 12233 | 15681 | 3.838e-03 | 3.054e-02 | 2.735e-04 | 5.624e-03 | 3.838e-03 | 3.055e-02 | 1.125e-04 | 1.773e-03 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | 10 | 14275 | 12233 | 15681 | 1.347e-04 | 4.141e-03 | 6.353e-05 | 1.977e-03 | 1.404e-04 | 4.409e-03 | 2.971e-05 | 7.993e-04 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | 200 | 14275 | 12233 | 15681 | 4.145e-03 | 2.758e-02 | 3.489e-04 | 9.337e-03 | 4.144e-03 | 2.756e-02 | 1.334e-04 | 2.946e-03 |
| SA 0.01-20 ConsGalerkinProj Drift Seed2 | 395 | 14275 | 12233 | 15681 | 3.995e-03 | 2.870e-02 | 2.431e-04 | 5.743e-03 | 3.995e-03 | 2.868e-02 | 1.050e-04 | 2.337e-03 |
| SA 0.01-20 NN Drift | 10 | 14275 | 12233 | 15681 | 1.052e-04 | 4.404e-03 | 6.725e-06 | 2.077e-04 | 1.123e-04 | 4.688e-03 | 2.650e-06 | 7.510e-05 |
| SA 0.01-20 NN Drift | 200 | 14275 | 12233 | 15681 | 4.459e-03 | 3.266e-02 | 3.694e-04 | 8.924e-03 | 4.459e-03 | 3.262e-02 | 1.724e-04 | 4.839e-03 |
| SA 0.01-20 NN Drift | 395 | 14275 | 12233 | 15681 | 3.838e-03 | 3.012e-02 | 3.265e-04 | 6.143e-03 | 3.838e-03 | 3.014e-02 | 1.310e-04 | 2.343e-03 |
| SA 0.01-20 NN Drift Seed2 | 10 | 14275 | 12233 | 15681 | 1.149e-04 | 4.188e-03 | 7.555e-06 | 1.947e-04 | 1.227e-04 | 4.466e-03 | 3.621e-06 | 1.063e-04 |
| SA 0.01-20 NN Drift Seed2 | 200 | 14275 | 12233 | 15681 | 4.159e-03 | 2.755e-02 | 1.934e-04 | 5.472e-03 | 4.158e-03 | 2.753e-02 | 8.664e-05 | 2.077e-03 |
| SA 0.01-20 NN Drift Seed2 | 395 | 14275 | 12233 | 15681 | 3.980e-03 | 2.845e-02 | 1.580e-04 | 4.239e-03 | 3.981e-03 | 2.843e-02 | 7.551e-05 | 2.307e-03 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 200 | 27780 | 12233 | 2176 | 8.860e-08 | 2.394e-06 | 8.052e-05 | 5.647e-03 | 6.195e-08 | 1.614e-06 | 8.285e-05 | 5.839e-03 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 250 | 27780 | 12233 | 2176 | 1.420e-06 | 2.591e-05 | 1.454e-04 | 8.677e-03 | 1.432e-06 | 2.591e-05 | 1.446e-04 | 8.679e-03 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 300 | 27780 | 12233 | 2176 | 2.261e-06 | 2.836e-05 | 1.455e-04 | 8.677e-03 | 2.288e-06 | 2.852e-05 | 1.447e-04 | 8.679e-03 |
| SA 0.05-0.10 ConsGalerkinProj OneShot | 395 | 27780 | 12233 | 2176 | 2.644e-05 | 1.931e-04 | 1.498e-04 | 8.698e-03 | 2.652e-05 | 1.935e-04 | 1.491e-04 | 8.700e-03 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 200 | 28819 | 12233 | 1137 | 8.538e-09 | 2.227e-07 | 1.228e-04 | 1.018e-02 | 8.552e-09 | 2.226e-07 | 1.266e-04 | 1.040e-02 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 250 | 28819 | 12233 | 1137 | 2.251e-06 | 1.486e-05 | 3.206e-04 | 1.736e-02 | 2.294e-06 | 1.497e-05 | 3.189e-04 | 1.737e-02 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 300 | 28819 | 12233 | 1137 | 4.768e-06 | 4.200e-05 | 3.207e-04 | 1.736e-02 | 4.820e-06 | 4.216e-05 | 3.190e-04 | 1.737e-02 |
| SA 0.10-0.15 ConsGalerkinProj OneShot SF2 | 395 | 28819 | 12233 | 1137 | 2.702e-05 | 1.697e-04 | 3.215e-04 | 1.737e-02 | 2.707e-05 | 1.700e-04 | 3.198e-04 | 1.738e-02 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 200 | 28819 | 12233 | 1137 | 4.636e-09 | 1.113e-07 | 6.109e-05 | 5.011e-03 | 4.201e-09 | 1.113e-07 | 6.293e-05 | 5.118e-03 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 250 | 28819 | 12233 | 1137 | 7.679e-07 | 8.962e-06 | 1.566e-04 | 8.466e-03 | 7.791e-07 | 8.961e-06 | 1.555e-04 | 8.471e-03 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 300 | 28819 | 12233 | 1137 | 1.553e-06 | 1.439e-05 | 1.566e-04 | 8.467e-03 | 1.575e-06 | 1.444e-05 | 1.555e-04 | 8.472e-03 |
| SA 0.10-0.15 ConsGalerkinProj OneShot | 395 | 28819 | 12233 | 1137 | 8.968e-06 | 5.544e-05 | 1.567e-04 | 8.465e-03 | 8.981e-06 | 5.554e-05 | 1.556e-04 | 8.470e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 400 | 29432 | 12233 | 524 | 1.915e-09 | 5.566e-08 | 2.856e-05 | 1.980e-03 | 1.711e-09 | 5.565e-08 | 2.907e-05 | 2.027e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 450 | 29432 | 12233 | 524 | 4.478e-07 | 3.785e-06 | 7.392e-05 | 3.882e-03 | 4.514e-07 | 3.785e-06 | 7.336e-05 | 3.878e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 500 | 29432 | 12233 | 524 | 8.623e-07 | 6.234e-06 | 7.392e-05 | 3.882e-03 | 8.647e-07 | 6.233e-06 | 7.336e-05 | 3.878e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot Long | 795 | 29432 | 12233 | 524 | 3.032e-06 | 1.764e-05 | 7.396e-05 | 3.881e-03 | 3.031e-06 | 1.764e-05 | 7.339e-05 | 3.877e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 200 | 29432 | 12233 | 524 | 1.468e-09 | 5.567e-08 | 2.379e-05 | 1.590e-03 | 1.453e-09 | 5.566e-08 | 2.388e-05 | 1.612e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 250 | 29432 | 12233 | 524 | 3.886e-07 | 5.511e-06 | 5.181e-05 | 3.236e-03 | 3.910e-07 | 5.510e-06 | 5.099e-05 | 3.221e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 300 | 29432 | 12233 | 524 | 7.409e-07 | 6.290e-06 | 5.182e-05 | 3.235e-03 | 7.435e-07 | 6.289e-06 | 5.100e-05 | 3.219e-03 |
| SA 0.20-0.25 ConsGalerkinProj OneShot | 395 | 29432 | 12233 | 524 | 2.266e-06 | 1.308e-05 | 5.184e-05 | 3.238e-03 | 2.265e-06 | 1.308e-05 | 5.102e-05 | 3.222e-03 |

**SST: node-wise discrepancy on never-displaced nodes** (`rans_turb_state_norms_SST.csv`). Nodes whose coordinates coincide with the reference to 1e-12 are compared directly; "band" = wall distance < 0.01 (the wall-layer nodes), "outer" = the other fixed nodes. Values are RMS and max of |case − reference| divided by the reference maximum over the same node set.

| case | step | fixed nodes | of which band | moved nodes | Eddy_Viscosity band RMS | band max | outer RMS | outer max | Turb_Kin_Energy band RMS | band max | outer RMS | outer max | Omega band RMS | band max | outer RMS | outer max |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SST 0.01-20 ConsGalerkinProj Drift | 10 | 14275 | 12233 | 15681 | 1.068e-03 | 2.211e-02 | 5.047e-05 | 1.386e-03 | 5.899e-05 | 2.281e-03 | 1.715e-06 | 5.956e-05 | 4.250e-07 | 9.666e-06 | 5.936e-08 | 1.464e-06 |
| SST 0.01-20 ConsGalerkinProj Drift | 200 | 14275 | 12233 | 15681 | 4.230e-03 | 8.556e-02 | 3.746e-04 | 1.199e-02 | 1.193e-03 | 9.996e-03 | 1.111e-05 | 3.460e-04 | 2.031e-06 | 1.047e-04 | 3.333e-07 | 9.529e-06 |
| SST 0.01-20 ConsGalerkinProj Drift | 395 | 14275 | 12233 | 15681 | 4.328e-03 | 6.440e-02 | 2.452e-04 | 4.854e-03 | 1.060e-03 | 1.078e-02 | 7.193e-06 | 1.094e-04 | 2.044e-06 | 5.362e-05 | 2.312e-07 | 3.082e-06 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | 10 | 14275 | 12233 | 15681 | 1.133e-03 | 2.755e-02 | 7.100e-05 | 1.450e-03 | 6.463e-05 | 2.322e-03 | 2.033e-06 | 5.205e-05 | 4.208e-07 | 1.197e-05 | 7.213e-08 | 1.478e-06 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | 200 | 14275 | 12233 | 15681 | 3.915e-03 | 7.338e-02 | 2.756e-04 | 4.859e-03 | 1.176e-03 | 1.059e-02 | 9.134e-06 | 1.920e-04 | 1.461e-06 | 5.592e-05 | 2.787e-07 | 4.978e-06 |
| SST 0.01-20 ConsGalerkinProj Drift Seed2 | 395 | 14275 | 12233 | 15681 | 4.382e-03 | 1.002e-01 | 2.431e-04 | 5.769e-03 | 1.165e-03 | 1.251e-02 | 7.078e-06 | 1.701e-04 | 2.267e-06 | 1.057e-04 | 2.227e-07 | 5.163e-06 |
| SST 0.01-20 NN Drift | 10 | 14275 | 12233 | 15681 | 1.046e-03 | 2.262e-02 | 5.353e-06 | 1.746e-04 | 5.177e-05 | 2.354e-03 | 1.827e-07 | 5.163e-06 | 4.202e-07 | 9.896e-06 | 6.841e-09 | 1.724e-07 |
| SST 0.01-20 NN Drift | 200 | 14275 | 12233 | 15681 | 4.197e-03 | 8.394e-02 | 3.934e-04 | 1.143e-02 | 1.224e-03 | 1.036e-02 | 1.177e-05 | 3.281e-04 | 2.047e-06 | 1.038e-04 | 3.428e-07 | 8.851e-06 |
| SST 0.01-20 NN Drift | 395 | 14275 | 12233 | 15681 | 4.295e-03 | 6.380e-02 | 2.746e-04 | 4.363e-03 | 1.086e-03 | 1.067e-02 | 8.990e-06 | 1.709e-04 | 2.015e-06 | 5.144e-05 | 2.842e-07 | 5.138e-06 |
| SST 0.01-20 NN Drift Seed2 | 10 | 14275 | 12233 | 15681 | 1.105e-03 | 2.687e-02 | 8.453e-06 | 2.601e-04 | 5.232e-05 | 2.059e-03 | 2.482e-07 | 7.552e-06 | 4.135e-07 | 1.197e-05 | 8.469e-09 | 2.179e-07 |
| SST 0.01-20 NN Drift Seed2 | 200 | 14275 | 12233 | 15681 | 3.859e-03 | 7.026e-02 | 2.009e-04 | 4.638e-03 | 1.208e-03 | 1.055e-02 | 5.925e-06 | 1.369e-04 | 1.452e-06 | 5.454e-05 | 1.871e-07 | 3.559e-06 |
| SST 0.01-20 NN Drift Seed2 | 395 | 14275 | 12233 | 15681 | 4.314e-03 | 1.007e-01 | 1.820e-04 | 5.901e-03 | 1.208e-03 | 1.248e-02 | 5.350e-06 | 1.738e-04 | 2.289e-06 | 1.060e-04 | 1.763e-07 | 5.243e-06 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 200 | 27780 | 12233 | 2176 | 6.321e-07 | 1.925e-05 | 8.470e-05 | 3.459e-03 | 5.999e-08 | 1.793e-06 | 6.741e-05 | 4.431e-03 | 6.239e-09 | 2.301e-07 | 4.317e-05 | 3.943e-03 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 250 | 27780 | 12233 | 2176 | 3.241e-06 | 4.774e-05 | 2.425e-04 | 8.246e-03 | 4.565e-06 | 5.787e-05 | 1.090e-04 | 5.972e-03 | 3.240e-08 | 9.205e-07 | 7.198e-05 | 5.292e-03 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 300 | 27780 | 12233 | 2176 | 3.453e-06 | 4.482e-05 | 2.461e-04 | 8.250e-03 | 1.040e-05 | 1.522e-04 | 1.090e-04 | 5.970e-03 | 4.179e-08 | 8.055e-07 | 7.203e-05 | 5.290e-03 |
| SST 0.05-0.10 ConsGalerkinProj OneShot | 395 | 27780 | 12233 | 2176 | 9.636e-06 | 8.875e-05 | 2.500e-04 | 8.292e-03 | 1.391e-05 | 9.106e-05 | 1.198e-04 | 5.975e-03 | 6.030e-08 | 1.227e-06 | 7.442e-05 | 5.295e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 200 | 28819 | 12233 | 1137 | 4.217e-08 | 1.252e-06 | 1.687e-04 | 9.471e-03 | 7.182e-09 | 1.015e-07 | 8.413e-05 | 7.755e-03 | 2.087e-09 | 1.151e-07 | 4.710e-05 | 3.431e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 250 | 28819 | 12233 | 1137 | 2.472e-06 | 2.664e-05 | 6.627e-04 | 2.926e-02 | 3.346e-06 | 3.472e-05 | 1.779e-04 | 9.824e-03 | 4.536e-08 | 1.381e-06 | 1.005e-04 | 5.887e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 300 | 28819 | 12233 | 1137 | 4.020e-06 | 3.990e-05 | 6.625e-04 | 2.924e-02 | 1.146e-05 | 1.270e-04 | 1.780e-04 | 9.827e-03 | 5.720e-08 | 1.496e-06 | 1.006e-04 | 5.887e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot SF2 | 395 | 28819 | 12233 | 1137 | 1.243e-05 | 8.976e-05 | 6.628e-04 | 2.921e-02 | 2.298e-05 | 9.991e-05 | 1.798e-04 | 9.825e-03 | 8.593e-08 | 2.305e-06 | 1.013e-04 | 5.889e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 200 | 28819 | 12233 | 1137 | 1.377e-08 | 3.576e-07 | 8.816e-05 | 4.764e-03 | 3.376e-09 | 6.768e-08 | 4.144e-05 | 3.815e-03 | 1.805e-09 | 1.151e-07 | 2.311e-05 | 1.682e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 250 | 28819 | 12233 | 1137 | 7.321e-07 | 1.019e-05 | 3.855e-04 | 1.756e-02 | 1.457e-06 | 1.790e-05 | 8.661e-05 | 4.874e-03 | 1.687e-08 | 4.603e-07 | 4.830e-05 | 2.770e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 300 | 28819 | 12233 | 1137 | 1.341e-06 | 1.445e-05 | 3.857e-04 | 1.755e-02 | 3.213e-06 | 3.804e-05 | 8.662e-05 | 4.872e-03 | 2.065e-08 | 4.603e-07 | 4.832e-05 | 2.770e-03 |
| SST 0.10-0.15 ConsGalerkinProj OneShot | 395 | 28819 | 12233 | 1137 | 4.516e-06 | 3.263e-05 | 3.861e-04 | 1.754e-02 | 7.888e-06 | 4.948e-05 | 8.679e-05 | 4.857e-03 | 3.074e-08 | 8.194e-07 | 4.842e-05 | 2.771e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 400 | 29432 | 12233 | 524 | 3.353e-09 | 1.192e-07 | 7.159e-05 | 5.175e-03 | 1.564e-09 | 6.769e-08 | 1.529e-05 | 1.044e-03 | 1.040e-09 | 1.151e-07 | 9.850e-06 | 5.672e-04 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 450 | 29432 | 12233 | 524 | 4.041e-07 | 4.351e-06 | 1.974e-04 | 1.437e-02 | 5.446e-07 | 6.430e-06 | 3.646e-05 | 1.807e-03 | 1.153e-08 | 3.452e-07 | 2.234e-05 | 1.005e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 500 | 29432 | 12233 | 524 | 5.997e-07 | 4.590e-06 | 1.975e-04 | 1.437e-02 | 2.145e-06 | 2.423e-05 | 3.647e-05 | 1.807e-03 | 1.208e-08 | 3.452e-07 | 2.235e-05 | 1.005e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot Long | 795 | 29432 | 12233 | 524 | 2.385e-06 | 1.168e-05 | 1.977e-04 | 1.438e-02 | 2.815e-06 | 3.581e-05 | 3.657e-05 | 1.807e-03 | 1.693e-08 | 3.452e-07 | 2.238e-05 | 1.005e-03 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 200 | 29432 | 12233 | 524 | 3.645e-09 | 5.960e-08 | 6.502e-05 | 4.387e-03 | 1.316e-09 | 6.768e-08 | 1.233e-05 | 8.941e-04 | 6.905e-11 | 7.192e-09 | 7.780e-06 | 4.530e-04 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 250 | 29432 | 12233 | 524 | 3.204e-07 | 5.424e-06 | 1.473e-04 | 8.588e-03 | 4.137e-07 | 4.399e-06 | 2.567e-05 | 1.477e-03 | 8.233e-09 | 1.151e-07 | 1.609e-05 | 8.867e-04 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 300 | 29432 | 12233 | 524 | 6.871e-07 | 7.152e-06 | 1.473e-04 | 8.588e-03 | 1.510e-06 | 1.638e-05 | 2.570e-05 | 1.476e-03 | 8.681e-09 | 1.151e-07 | 1.610e-05 | 8.862e-04 |
| SST 0.20-0.25 ConsGalerkinProj OneShot | 395 | 29432 | 12233 | 524 | 1.126e-06 | 7.000e-06 | 1.474e-04 | 8.588e-03 | 1.939e-06 | 1.347e-05 | 2.580e-05 | 1.477e-03 | 9.513e-09 | 1.151e-07 | 1.613e-05 | 8.868e-04 |

## III.7 Figures and files

Figures (`Figures/RANS/`):

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

Tables (`RANS_PostProcessing/`):

- `rans_deformed_reference.csv`
- `rans_deformed_surface.csv`
- `rans_drift_series.csv`
- `rans_drift_trends.csv`
- `rans_oneshot_fit.csv`
- `rans_oneshot_lags.csv`
- `rans_oneshot_series.csv`
- `rans_pressure_probe.csv`
- `rans_surface_forces.csv`
- `rans_surface_norms.csv`
- `rans_surface_profiles.csv`
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
