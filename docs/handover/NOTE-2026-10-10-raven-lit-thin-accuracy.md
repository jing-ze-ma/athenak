# Thin-region accuracy vs. thick-limit AP in implicit VET / DOM schemes: what the papers actually do

Literature study for the AthenaK `<rad_m1>` implicit VET face flux (xthinfix-1010). Written 2026-10-10.

## 0. Sources and how they were read

I read the full arXiv LaTeX sources (downloaded with `arxiv.org/e-print/<id>`) of:

| short | paper | arXiv | what was read |
| --- | --- | --- | --- |
| JSD12 | Jiang, Stone & Davis 2012, ApJS 199, 14 | https://arxiv.org/abs/1201.2223 | Sect. 3 (algorithm), 3.3 (radiation subsystem), 4-6 (tests), 9 (summary), App. A (the matrix) |
| JSD13 | Jiang, Stone & Davis 2013, ApJ 767, 148 ("Saturation of the MRI in strongly radiation-dominated disks") | https://arxiv.org/abs/1303.1823 | **App. A, "Improvements to the numerical algorithm"**: the tau-dependent HLLE speed for the implicit VET moment solver. This is the origin of the VETTAM correction. Jiang 2021 cites it as "Jiangetal2013b" |
| JSD14 | Jiang, Stone & Davis 2014, ApJS 213, 7 (explicit time-dependent intensity solver) | https://arxiv.org/abs/1403.6126 | Sect. 4.2 (transport step) and the beam, dynamic diffusion and linear wave tests. It is explicit, and short characteristics are not used here (that is Davis, Stone & Jiang 2012). It is relevant only because it gives the alpha(tau) transport speed |
| J21 | Jiang 2021, ApJS 253, 49 (implicit discrete ordinates) | https://arxiv.org/abs/2102.02212 | Sect. 3 complete, Sect. 4 complete, Sect. 5, **App. A complete** |
| VETTAM | Menon et al. 2022, MNRAS 512, 401 | https://arxiv.org/abs/2202.08778 | Sect. 2.2 complete, Sect. 3.1-3.6, Sect. 4.3 (caveats), **App. A complete** |

**Equation numbers** were computed from the arXiv LaTeX with a numbering script (`\nonumber` and cases-aware, appendices restart at A1). I spot-checked them against ar5iv for J21: eq. (17) = `f = 1 - exp(-tau_c^2)` and eq. (19) = the flux divergence, both correct. The journal versions can differ by one in places. Section numbers come from the source.

**Not verified from primary text (flagged):**
- Sekora & Stone 2010 (SS10), eq. 39 (the HLLE flux all four use) and its Sect. 4.2 (oscillation and monotonicity argument against higher-order implicit schemes). I inferred the HLLE form of JSD12 from its App. A coefficients (Sect. 1.1).
- Audit et al. 2002.
- The figures. All orders and errors quoted below are the papers' own text statements. I did not read slopes off the plots.
- AREPO-IDORT (Ma et al. 2025). Our `idort` weight is attributed to it only by our code comment. I did not re-verify it here.

---

## 1. Paper by paper

### 1.1 JSD12: the original implicit VET moment scheme (Athena)

1. **Spatial discretisation.**
   - The radiation subsystem eq. (8) (`dE/dt + C div F = C S_E`, `dF/dt + C div(f E) = C S_F`) uses **first-order (dc) reconstruction**: "Since the backward Euler differencing is only first-order in time, first-order spatial reconstruction is used to compute the left- and right-states for the HLLE solver" (Sect. 3.3.2, after eqs. 20-21).
   - It reconstructs the conserved `E_r, F_r`, piecewise constant. There is no limiter and no lagged high-order part. The gas uses 2nd/3rd-order primitive reconstruction, which is separate.
2. **Face flux.**
   - Plain HLLE (SS10 eq. 39), fully implicit: eqs. (20)-(21), with the 3-D matrix in App. A (`theta_0..15`, `phi`, `psi`, `varphi`).
   - From the coefficients (`Ci0 = (sqrt f_i - sqrt f_{i-1})/(sqrt f_i + sqrt f_{i-1})`, `theta_4 = -d1x (1+Ci0) sqrt f_{i-1}`, `theta_5 = -d1x (1+Ci0)`, ...), the face speeds are `S_L = -C sqrt(f_xx(i-1))` and `S_R = +C sqrt(f_xx(i))`, i.e. the cell-centred VET diagonal of each neighbour. I checked that `F_E = (S_R F_L - S_L F_R + S_L S_R (E_R - E_L))/(S_R - S_L)` reproduces `theta_4` and `theta_5`.
   - **There is no optical-depth modification in JSD12.**
3. **Thin vs thick.** No switch. At full speed `+-C`, a beam with `f_xx = 1` gets the exact upwind flux (`F = cE`). In thick cells HLLE adds `D_num ~ C dx sqrt(f)/2`. JSD13 later recognised this as wrong for scattering-dominated thick media; it is invisible in JSD12 because the absorption source dominates there.
4. **Time.**
   - Radiation: backward Euler. "higher-order implicit time integration schemes can lead to oscillatory solutions with large time steps. Thus, to ensure a non-oscillatory method, we restrict the update to first-order backward Euler" (Sect. 3, citing SS10).
   - Gas: explicit, 2nd-order modified Godunov. The split order is arbitrary.
   - The VET comes from short characteristics at step start (Step 1) and is held fixed. The source terms are made linear (T estimated, gas held fixed), so there is **one linear solve per step**.
   - Solvers (Sect. 8.2): LIS/Hypre GMRES with BoomerAMG, and their own multigrid, which is the "workhorse".
5. **Tests and orders.**
   - Linear waves (Sect. 6.1): "the radiation dominated case converges at first order, while the gas pressure dominated case converges at second order", for both RHD and RMHD. Damping and phase are inaccurate for `tau_a > C/P` (Sect. 9).
   - Marshak, tophat (crooked pipe), shadow (two beams at +-14 deg, umbra and penumbra): qualitative only. The authors argue the first-order BE+VET still "can maintain very sharp gradients".
   - **No spatial convergence study in an optically thin region.**
6. **Failure modes mentioned.**
   - Higher-order implicit time integration oscillates.
   - Energy non-conservation (`R = 0.05` fix).
   - The time-independent VET makes the tophat heat too early.
   - Inaccurate linear waves at `tau_a > C/P`.
   - Nothing about negative E or superluminal F.

### 1.2 JSD13 App. A: the tau-dependent HLLE speed (the VET "AP fix")

1. **Spatial discretisation.** Still dc HLLE for `(E, F)` in the implicit BE matrix. **New:** the advective piece `div(v E_r)` is split off and done **explicitly with a 2nd-order van Leer slope** of `E_r^n` (App. A.2; face velocity is the mean of the two cells). The `div(v.P_r)` piece is implicit and centred (App. A.3, eq. A7). This is the first "high-order part lagged/explicit, low-order part implicit" construction in this line of work.
2. **Face flux.**
   - The HLLE max/min speeds are replaced by `+-C*_eff`:
     ```
     C_eff    = sqrt(f) C sqrt(1 - sigma_t^2/(4 f k^2))                        (A2)  [dispersion of the telegraph system A1]
     C*_eff   = C sqrt( (f/tau) [1 - exp(-tau)] ),  tau = (10 dl sigma_t)^2 / (2 f)    (A3)
     ```
   - `sigma_t = sigma_aF + sigma_sF` (flux-mean absorption plus scattering), `dl` = cell size, `f` = the VET component along the face normal. **On the face, `sigma_t` is the arithmetic mean of the two cells**, chosen "to make the implicit matrix easier to invert for the cases with sharp opacity jump". The same speed enters the E and F rows (the whole HLLE vector).
   - Wavelength: the speed is evaluated at `k = 1/(10 dl)` ("a wave length of ten cells").
   - **They state explicitly that the physically "right" `k = 1/(dt C)` was rejected because it "makes the resulting matrix from the implicit radiation subsystem very hard to invert".** That is a dt-dependent switch rejected for conditioning, the same failure axis as our `xthin` override.
3. **Thin vs thick.**
   - Smooth, opacity-only, resolution-dependent: `C* -> sqrt(f) C (1 - tau/4)` for `tau -> 0`, and `C* -> C f sqrt2/(10 tau_cell)` for `tau_cell >> 1`.
   - The transition (`tau = 1`) sits at `tau_cell = sqrt(2f)/10`: 0.08 for f = 1/3, 0.14 for f = 1.
   - `C*/C sqrt f` = 0.72 at `tau_cell = 0.1`, 0.27 at 0.3 and 0.08 at 1 (f = 1/3). So the HLLE is already 30 % "less upwind" at `tau_cell = 0.1`.
   - There is no gradient or Knudsen criterion and no threshold.
4. **Time.** BE as in JSD12 plus the explicit van Leer advection of `E` (CFL-limited by the flow speed).
5. **Tests.** Dynamic diffusion with pure scattering at `sigma_t = 4e4, 40, 1`, `C = 514.4`, `v = 1` (App. A.4, eqs. A9-A10): "agree ... very well". Without the fix, "the diffusion time is too short". No orders are given.
6. **Failure modes.**
   - Plain HLLE numerical diffusion dominates in scattering-dominated thick media.
   - First-order advection is too diffusive.
   - The dt-based wavenumber gives an ill-conditioned matrix.

### 1.3 JSD14 (explicit, for the alpha(tau) only)

- Transport speed `alpha C` with `alpha = sqrt((1 - exp(-tau))/tau)`, `tau = [10 dl (sigma_a + sigma_s)]^2` (eq. 14, Sect. 4.2), applied to `C n.grad I~`.
- The `3 n.v J` advective part is split off (eq. 13) and moved at speed `|v|`.
- Both parts use the explicit **2nd-order van Leer** upwind interpolation (Stone & Mihalas 1992, eq. 4).
- Linear waves: 2nd order at `tau_a = 0.01`, 1st order at `tau_a = 10`, `P = 10` (Sect. 5, "Linear Wave Convergence"; attributed to operator splitting).
- Crossing beams: "a small amount of numerical diffusion ... reduced with increasing resolution". No order is given.
- Without the alpha fix, numerical diffusion dominates at `tau_s = 625`.

### 1.4 J21: implicit discrete ordinates (Athena++)

1. **Spatial discretisation.**
   - Intensities `I_n` per angle, **first-order (dc) left/right states** `I(i-1), I(i)` (Sect. 3.2.1).
   - Stated reasons (Sect. 3.2): "higher order spatial reconstruction without flux limiter cannot guarantee total variation diminishing ... The commonly used flux limiters will make the transport term non-linear, which is very hard to solve implicitly." And again: "In principle, these numerical errors can be minimized with second or third order spatial reconstructions. However, it is not an option for our implicit scheme."
   - **The only 2nd-order piece is lagged/explicit.** The transport is split as `c div(n I) = div[(c n - f v) I] + div(f v I)` (eq. 16), with `f = 1 - exp(-tau_c^2)` (eq. 17). The `div(f v^m I^m)` term uses step-start intensities, "second order spatial reconstruction for `I_n^m`" and an upwind flux on `f v^m` (eq. 18, App. A eq. A1). It is CFL-safe because `dt` satisfies the flow CFL.
2. **Face flux (App. A).**
   - Modified HLLE per angle, with the advective velocity `a = c mu_x - f v_x^m(i-1/2)` (eq. A3):
     ```
     F_n(i-1/2) = S+/(S+ - S-) a I_L - S-/(S+ - S-) a I_R + S+ S-/(S+ - S-) (I_R - I_L)          (A3)
     mu_x > 0:  S+ =  c mu_x sqrt[(1 - e^{-tau_c^2})/tau_c^2],  S- = -c mu_x sqrt[(1 - e^{-tau_c^4})/tau_c^2]   (A4, A5)
     mu_x < 0:  mirror image
     tau_c = alpha [rho(i-1) + rho(i)] [kappa_a(i-1) + kappa_a(i) + kappa_s(i-1) + kappa_s(i)] dx,  alpha = 5 (default)
     ```
   - Uniform medium: `tau_c = 20 tau_cell`.
   - Limits: `tau_c -> 0` gives the pure upwind flux (`S- ~ -c mu tau_c -> 0`). `tau_c -> inf` gives `F -> (a/2)(I_L + I_R) - (c|mu|/(2 tau_c))(I_R - I_L)` (eq. A6).
   - "alpha is chosen to minimize numerical diffusion while still keep the numerical scheme stable ... the solution is independent of alpha as long as it is large enough."
   - Written as a convex blend (my algebra, with `p = S+/(c mu)`, `q = -S-/(c mu)`): `F = F_central + w_J (F_upwind - F_central)`, with `w_J = (p - q + 2pq)/(p + q)`. So `w_J -> 1` in thin cells and `w_J -> 1/tau_c = 1/(20 tau_cell)` in thick cells. Values: `w_J` = 0.80 at `tau_cell = 0.05`, 0.49 at 0.1, 0.17 at 0.3 and 0.05 at 1.
   - The `tau^4` on the downwind speed makes the upwind-to-central transition smooth (`S-/S+ ~ tau_c` for small `tau_c`).
3. **Thin vs thick.** Smooth, opacity-only (`tau_c` per face from the two cells' rho and kappa), resolution-dependent. No threshold, no gradient criterion. The advective split weight `f = 1 - exp(-tau_c^2)` uses the same `tau_c`.
4. **Time and solver.**
   - BE for all transport and source terms, with `T^{m+1}` nonlinear (eq. 14). Opacities, `rho` and `v` are held at step `m`.
   - The "Jacobi-like" iteration (eq. 21) does a local per-cell solve over all angles, reduced to a quartic in `T`, with neighbours lagged by one iteration. It stops at `Delta I < 1e-5` (eq. 26).
   - Robustness (App. A, last paragraph): convergence "depends strongly on the values of g1". The alternative splits `g1 = g1' (positive) + g1'' (negative)` and lags `g1''`. It is "much more robust", but slower when both converge.
   - Coupling to the VL2 hydro: radiation solve over `dt/2` in the predictor, then over `dt` in the corrector, each BE (Sect. 3.3).
   - There is no VET (DOM). J21 notes that VET schemes effectively use the Eddington tensor at the start of the step, while DOM uses the end-of-step one (Sect. 4.3).
5. **Tests and orders.**
   - Linear waves (Sect. 4.7): "When the error is dominated by the RT module, it only shows first order convergence". Part of that truncation error comes from "different upwind directions for specific intensities propagating along opposite directions". At `N = 512`, "increasing tau_c ... by a factor of 2 from our default value can help decrease the L1 error" (a symptom of a residual, resolution-independent numerical diffusivity).
   - Crossing beams in vacuum (Sect. 4.2): "very similar to the case with first order spatial reconstruction but more diffusive compared with ... second order" (of JSD14). It needed 100 iterations in step 1 to reach 1e-3.
   - Static and dynamic diffusion (Sect. 4.4): match the analytic solution. "we need to use a larger value of the parameter tau_c to get the same accuracy with a larger kappa_s".
   - Non-LTE atmosphere (Sect. 4.5): matches. "The relative error is also larger in the optically thin surface and the error increases with reducing epsilon for a fixed spatial resolution". No order is given.
   - Homogeneous sphere (Sect. 4.6): sharp `r = 1` transition captured "without any artificial numerical effect".
   - **No thin-region spatial convergence ladder.**
6. **Failure modes.**
   - Upwind fluxes overwhelm the physical diffusion in thick cells (the two-angle argument in Sect. 3.2.1).
   - The iteration can converge slowly or diverge when `g1` is not diagonally dominant.
   - Vacuum transport is the worst case for the iteration.
   - First-order beams are diffusive.
   - Ray effects.
   - Nothing about negative intensities or superluminal flux.

### 1.5 VETTAM (Menon et al. 2022, FLASH, AMR)

1. **Spatial discretisation.** "piecewise constant (first-order) reconstruction, using the state of the conserved quantities `U_r` at time `t^{n+1}`" (Sect. 2.2.3). This is inside the implicit operator. Nothing is lagged at higher order.
2. **Face flux.**
   - HLLE (SS10 eq. 39) with
     ```
     C_HLLE = sqrt(f) c sqrt[(1 - e^{-tau_c})/tau_c]          (34)
     tau_c  = (10 dl rho kappa_R)^2 / (2 f)                    (35)
     tau_c at the face = (tau_L + tau_R)/2 on a same-level face; upstream tau (tau_L for the right-going wave, tau_R for the left-going one) at an AMR level jump   (36)
     ```
   - This is JSD13 eq. (A3) verbatim. One difference: VETTAM averages `tau_c`, which is proportional to sigma squared, and states that JSD13 does the same. JSD13's text actually says the **opacity** is averaged. The two differ at opacity jumps.
   - Coarse-fine conservation: the coarse flux is replaced by the mean of the fine fluxes inside the implicit system (eq. 37).
3. **Thin vs thick.** Same as JSD13 (smooth, opacity-only, transition `tau_cell ~ 0.08-0.14`).
4. **Time and solver.**
   - BE for `(E, E_r, F_r)`, explicitly "restrict[ed] ... to first-order because higher-order implicit time integration schemes have been found to lead to oscillatory solutions when using large time steps" (Sect. 2.2.2).
   - Picard on `T_*` (eqs. 38-41): each pass does a linearised BE solve with PETSc GMRES and ASM, then a local Newton for `T`. Tolerances: `eps_R = eps_N = 1e-6`, `eps_P = 1e-3`.
   - The VET comes from the hybrid-characteristics ray tracer (HEALPix, randomly rotated) at step start, fixed through the step.
   - The `O(v/c)` terms are explicit by default, with an implicit runtime switch for dynamic diffusion.
5. **Tests.**
   - Radiating pulse, `dt = dx/c`, 1024 cells (Sect. 3.1). Static diffusion L1: 5 %, 3.7 %, 3.3 %. Equilibrium diffusion: 2.4 %, 3 %, 3.2 %. **Streaming: 3.7 %, 16 %, 27 %, growing with time**: "the pulse has diffused ... This is not surprising considering that our scheme is only first-order accurate in space and time".
   - Dynamic diffusion: L1 4.5 %, 5.6 %, 5.8 %.
   - Advecting pulse: max 0.7 %.
   - App. A1: without the correction, the weak-diffusion pulse needs 4096 cells (`tau_cell ~ 0.1`) to converge. With it, 256 cells (`tau_cell ~ 1.5`) are already close.
   - App. A2: AMR level boundary. No correction gives a discontinuity in `F`, the averaged tau gives "spurious oscillations", and the upstream tau is smooth and matches UG2048.
   - Shadows from point and diffuse sources versus M1 and Eddington closures (Sect. 4.1).
   - **No convergence orders reported anywhere.**
6. **Failure modes (Sect. 4.3 caveats).**
   - Streaming is diffusive.
   - "implicit methods perform poorly when trying to capture the propagation of individual wave modes".
   - The step-start VET is inadequate when the radiation field changes within a step.
   - Picard converges only linearly.
   - Energy conservation is only to the solver tolerance.
   - Higher-order implicit monotonicity is an open problem (citing SS10 Sect. 4.2).
   - Averaged tau oscillates at AMR jumps.
   - Nothing about negative E or `|F| > cE`.

---

## 2. Side-by-side

| | JSD12 | JSD13 App. A | J21 | VETTAM | **ours** (xthinfix) |
| --- | --- | --- | --- | --- | --- |
| unknowns | `E, F` collocated | `E, F` collocated | `I_n` | `E, F` collocated (+ gas E) | `E` cell, `F0` on faces, eliminated (scalar E operator) |
| reconstruction in the implicit operator | dc | dc | dc | dc | dc (`implicit_hr_recon = dc`), or plm half-range with a frozen van Leer limiter (`plm`, Fix B) |
| lagged/explicit high-order piece | none | van Leer advection of `E` (explicit) | 2nd-order reconstruction of `f v I^m` (explicit) | none | plm deferred correction of the enthalpy flux (`implicit_enthalpy = plm`); plm limiter frozen per step (`qs_mode = step`) |
| thick-limit flux | plain HLLE (wrong for scattering) | centred `(2 dx)` plus HLLE at `eps ~ 1/tau_cell` | per-angle central plus `c|mu|/(2 tau_c)` | as JSD13 | compact face-eliminated central (exact harmonic face D) |
| `D_num/D_phys`, thick | `~ tau_cell` | **0.071 (f = 1/3), constant** | **0.0375, constant** | 0.071 | `0.75 w tau_cell`: **-> 0** for `idort` (`w ~ tau0/(alpha tau)^2`). Exactly 0 where `w = 0` |
| thin/thick switch | none | opacity only, smooth: `tau = (10 tau_cell)^2/(2f)` | opacity only, smooth: `tau_c = 20 tau_cell`, `e^{-tau^2}` / `e^{-tau^4}` | opacity only, smooth, upstream at AMR | opacity × beam measure (`idort_f`), plus the X override (`all`/`steep`, dt-dependent) or a Knudsen override (`kn`) |
| upwind fraction at `tau_cell` = 0.1 / 0.3 / 1 / 3 | 1 | 0.72 / 0.27 / 0.08 / 0.03 | 0.49 / 0.17 / 0.05 / 0.017 | 0.72 / 0.27 / 0.08 / 0.03 | `idort` alone: 0.90 / 0.72 / 0.33 / 0.08 (then × smoothstep(\|H\|/J; 0.3, 0.6) for `idort_f`) |
| time | BE | BE | BE (×2 per VL2 step) | BE | BE, with SDIRK2/TR-BDF2 options |
| solver | one linear solve per step (GMRES/AMG, MG) | same | nonlinear Jacobi + local quartic | Picard + GMRES/ASM | Picard + BiCGStab/line, Anderson optional |
| reported thin-region order | none | none | "first order" (beam, qualitative) | none (streaming L1 grows to 27 %) | measured X ladders (TASKS: sq has the right order, but diverges on beams) |

(The upwind fraction is the coefficient of `F_upwind - F_central` in the convex form. For the moment-HLLE schemes it is `eps = C*/(C sqrt f)`.)

---

## 3. What this means for the user's question

**(a) No published implicit VET/DOM scheme is better than first order in space in ANY regime.** All three keep the implicit operator dc for linearity and an M-matrix-like structure. J21 says explicitly that limiters make an implicit transport operator intractable. The only 2nd-order spatial pieces (JSD13, J21) are **explicit, CFL-limited by the gas velocity, and built from step-start data**. None of them lags a high-order part of the *light-speed* transport. None reports a convergence order in a thin region; their beam and streaming tests are qualitative or single-resolution.

**(b) Their AP is approximate.** With the speed scaled by `eps ~ 1/tau_cell`, the HLLE numerical diffusivity is a **constant fraction of the physical one**: 3.75 % for J21's `alpha = 5`, 7.1 % for JSD13/VETTAM at f = 1/3. It does not shrink with `dx` until `tau_cell` falls below about 0.05-0.1, where the correction switches off and the plain HLLE `O(dx)` error takes over. J21's own remarks confirm this resolution-independent error: "increase `tau_c` to get the same accuracy with a larger `kappa_s`", and "doubling `tau_c` lowers L1 at N = 512". The thick diffusion operator is also the collocated `2 dx` centred stencil, kept stable only by that residual dissipation. **Our compact face-eliminated central flux is strictly better in the thick limit**: exact in the limit, second order, the exact harmonic face diffusivity at opacity jumps, and no odd-even mode.

**(c) Any weight built from `tau_cell` alone goes to 1 (pure dc upwind) as `dx -> 0` at fixed physics.** That holds for J21's `w_J`, the JSD13/VETTAM `eps` and our `idort`. Under refinement every region therefore ends up in the dc upwind branch, so X order -> 1 *everywhere*, including physically thick and intermediate layers once `tau_cell < ~0.1-0.6`. To get order 2 in thick and intermediate regions under refinement, the weight must depend on something resolution-independent: the field's Knudsen number `|grad E|/(chi E)`, or the formal-solution flux factor `|H|/J`. Alternatively, the upwind branch itself must be 2nd order (plm). **None of the three papers does either.** This part is our own territory.

**(d) Conditioning.** JSD13 rejected a dt-dependent wavenumber (`k = 1/(c dt)`) because it made the implicit matrix "very hard to invert". The published resolution is to keep the switch **dt-independent and always slightly dissipative**. JSD13, J21 and VETTAM never let the dissipation reach exactly 0. J21 additionally lags the negative part of `g1` for robustness. VETTAM needed the upstream tau at AMR jumps to avoid oscillations, an opacity-switch analogue of a sharp transition in the weight.

## 4. Mapping to our options

| our option | published counterpart | where we differ |
| --- | --- | --- |
| `thick_flux = scaled` (explicit M1, `rad_m1_closure.hpp` `M1_THICK_SCALED`): `eps = sqrt[(1 - e^{-tc^2})/tc^2]`, `tc = scaled_prefactor tau_face`, both speeds, both fluxes | JSD13 eq. (A3) / VETTAM eqs. (34)-(35) **exactly** if `scaled_prefactor = 10/sqrt(2f)` (12.2 for f = 1/3, 7.1 for f = 1). J21 corresponds to `20` but is per-angle with the `tau^4` downwind asymmetry | our prefactor is f-independent (default 20 = J21), so a VETTAM-like run needs about 12. VETTAM averages `tau` (∝ σ²); we average σ, as JSD13 does. VETTAM's AMR rule is the upstream tau at level jumps (we have no AMR here). `rad_m1_design.md` Sect. 3 lists VETTAM and JSD14 as "only through summaries" / "NOT read"; that can now be upgraded: JSD13 App. A is the primary source |
| `thick_flux = ap_hll` (Berthon/Bloch alpha) | not in these papers | compact AP; none of the five uses it |
| `implicit_flux = central` (face-eliminated staggered `F0`) | **no counterpart.** The closest is the thick limit of J21 eq. (A6), but ours is compact and adds no dissipation | better than published in thick/intermediate regions. In thin regions it is the BE staggered moment system, which none of the papers uses there |
| `implicit_flux_beam = halfrange` upwind part `c (r+_L E_L + r-_R E_R)`, `r± = H±/J` from the VET rays | **the angular integral of J21's `tau -> 0` per-angle upwind flux** (eq. A3 with `S- = 0`), with the intensities replaced by the lagged formal-solution half-range ratios times the current E | J21 evolves the intensities, so its half-range moments are consistent with the end of the step. Ours freeze the angular shape at the VET/ray state. In transients with `X = c dt/dx >> 1` (a front crossing many cells per step) the lag is O(1). This is the VET-lag caveat both J21 (Sect. 4.3) and VETTAM (Sect. 4.3) name |
| `implicit_flux = blend`, `(1 - w) F_central + w F_hr` | J21 has the same convex structure (`F = F_c + w_J (F_up - F_c)`); for JSD13/VETTAM the weight is `eps` | J21's thick limit keeps `w_J ~ 1/tau_c` (constant `D_num` fraction). Ours goes to exactly 0 |
| `implicit_blend = idort`: `w = 1/(1 + x + x^2/tau0)`, `x = alpha tau_f` | J21 `w_J` (`-> 1/(20 tau_cell)`) and JSD13/VETTAM `eps` (`-> 0.08/tau_cell`) are the counterparts | **ours is far more upwind in intermediate regions:** `w = 0.5` at `tau_f = 0.62`, against 0.1 (J21) and 0.15 (VETTAM). `D_num/D_phys` is 16-25 % for `tau_cell` 0.3-3 (J21 3.75 %, VETTAM 7 %) but falls to 0 deep inside (truly AP, J21/VETTAM are not) |
| `implicit_blend = idort_f` (× smoothstep of \|H\|/J, `flo/fhi = 0.3/0.6`) | **no counterpart.** The papers' switches are opacity-only | `|H|/J ≈ Kn/3` in a diffusive field, so this is effectively a Knudsen gate at Kn ~ 1-2. It comes from the VET rays, i.e. **frozen per step**, so it costs nothing in Picard. This is the cheapest resolution-independent gate available, and it satisfies (c) |
| `xthin` `all`/`steep` (`X = (c dt/dx)/(1 + c dt chi_f)`, steep `X^8`) | **none.** JSD13 explicitly *rejected* a dt-dependent criterion (conditioning) | ours forces upwind on transparent faces to restore conditioning. The papers instead never drop the dissipation to 0 anywhere and keep the switch dt-free. A dt-dependent weight also makes the spatial scheme depend on the time step, which interferes with time-convergence (T) ladders |
| `xthin` `kn` (`R = |E_R - E_L|/(tau_f max E)` of the Picard iterate) | none (no gradient criteria in the three) | physically the right variable (resolution-independent, (c)). But evaluated on the **iterate**, it makes the operator nonlinear through Picard; the papers keep their weights fixed per step |
| `implicit_hr_recon = plm` (+ `qs`, `qs_mode = step`, `hr_pos = bound`) | the "explicit/lagged 2nd-order part" of JSD13 App. A.2 / J21 eq. (18) | they apply it only to the **flow-speed** advection (Courant number ≤ 1). Ours applies it to the light-speed half-range transport. That is stable only where the field is quasi-steady over the step, which is what `qs` enforces. `qs_mode = step` (limiter from `E^n`, row fixed through Picard) is the closest analogue of their "operator built once per step" |
| BE + Picard; VET lagged | all four | J21 lags nothing angular (DOM). VETTAM and JSD12 compute the VET at step start, as we do |
| 2nd order in time (SDIRK2/TR-BDF2 options) | none: all restrict themselves to BE because higher-order implicit schemes oscillate at large dt (SS10, quoted in JSD12 and VETTAM) | their objection is to non-L-stable schemes (Crank-Nicolson type). L-stable SDIRK2/TR-BDF2 is not ruled out by anything in these papers. SS10 was not read, so this is not verified |

## 5. Most promising for "2nd order (space+time) in thick/intermediate, >= 1st order in space in thin, cheap"

1. **Keep our compact central flux as the thick/intermediate operator.** None of the published approaches is second order there. All of them carry a resolution-independent `D_num` of 4-7 % and the `2 dx` stencil.
2. **Gate the upwind/half-range weight on a resolution-independent quantity that is frozen per step.** The published lesson (JSD13, J21, VETTAM) is to build the face weights once per step so the system is linear, and keep them dt-independent.
   - The best candidate in our code is `idort_f`'s `|H|/J` from the formal solution, which is already frozen with the VET.
   - The alternative is the Knudsen ratio computed from `E^n` (not the Picard iterate).
   - Either way, w stays 0 wherever the field is diffusive, however thin the cell, so X order 2 survives refinement. Where `Kn >~ 1` (true transport) the dc half-range upwind gives order 1, which meets the thin criterion.
   - This replaces the X-based override (`all`/`steep`), which is the one ingredient the literature explicitly warns against (JSD13: dt-dependent switch, ill-conditioned matrix).
3. **For conditioning at mixed faces, borrow the published trick rather than X:**
   - keep a small floor on the dissipation (`w >= w_min ~ 1/(20 tau_cell)` capped, as J21 does implicitly), or
   - lag the non-M-matrix (positive off-diagonal) part as J21 does with `g1''`, which our plm/`kill`/`bound` machinery partly emulates.
   - At sharp opacity/level jumps, take the weight from the **upstream** cell (VETTAM eq. 36). Averaging oscillated.
4. **Second order in the thin/beam part, if needed:** only by the lagged route (plm limiter frozen per step from `E^n`: `qs_mode = step`). It must be restricted to faces where it is a contraction (quasi-steady or `X <~ 1`), consistent with JSD13/J21 applying their explicit 2nd-order pieces only at Courant number ≤ 1. Nothing in print makes an implicit light-speed upwind scheme second order at `X >> 1`. The acceptance criterion only asks for order ≥ 1 there, so dc half-range is sufficient and cheapest.
5. **Time:** keep an L-stable 2nd-order scheme (TR-BDF2/SDIRK2). The papers' "BE only" rule targets oscillation from A-stable non-L-stable schemes. They report no attempt with L-stable ones.

## 6. Caveats

- Effective upwind fractions, `D_num/D_phys` and the `w_J` formula are my algebra from the published formulas (isotropic `<|mu|> = 1/2`, half-range `r± = ±1/4`, uniform medium). They are not numbers the papers give.
- The JSD12 HLLE speed assignment (`S_L` from `f(i-1)`, `S_R` from `f(i)`) was reverse-engineered from its App. A coefficients, because SS10 eq. 39 was not read.
- VETTAM's statement that JSD13 averages `tau_c` contradicts JSD13's own text, which averages the opacity. I report both as written.
- Equation numbers come from the arXiv sources, spot-checked for J21 only.
