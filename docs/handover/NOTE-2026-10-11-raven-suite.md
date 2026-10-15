# NOTE suite (suite-1011): the literature radiation tests on the merged production candidate (prodcand-1011)

Work dir: /raven/ptmp/jinma/suite_1011. MPCDF Raven, CPU (`general`, 72 cores, OpenMP, one MeshBlock).

**Question.** Does the latest recommended setting pass the tests of the IDORT paper (Ma, Pakmor, Justham & de Mink, arXiv 2503.16627), and the extra ones of Jiang 2021 (J21, 2102.02212) and Jiang, Stone & Davis 2012 (JSD12, 1201.2223)?

**Short answer.**
- **REC-be** (recommended keys, backward Euler) passes every previously reproduced test. No regression against the per-branch results.
- **REC** (recommended keys, hesdirk2 + `time2_vet_sc=predict`) passes everything except two cases. Both are caused by `time2_vet_sc=predict`, and both pass with `time2_vet_sc=lag`:
  - **thermal equilibrium case 2** creates a spurious momentum: v_x = 1.71 in a uniform periodic box, E 1 % low;
  - **crossing beams at c dt/dx = 256** are FATAL: Picard diverges at cycle 2.
  - Both bugs are already on the per-branch binaries (42aa6bbc, 9e1f2686, mq3). They are not merge regressions; earlier tests ran BE only. The predfix-1011 agent is fixing them.
- **Formerly "not done" tests:**
  - **radiation shocks M = 1.2 / 2 / 3: PASS** against the IDORT transport reference, with L1 ≤ 0.14 %.
  - **homogeneous sphere** (as a Cartesian cylinder): PASS at 512^2, marginal at 256^2. The error is the formal solution's angular smearing.
  - **crooked pipe:** reproduced. It converges in dt and is self-consistent, but there is no quantitative reference in Cartesian geometry.
  - **multi-source bubbles: FAIL.** The fronts lag the light cone, and total E is 0.4 of exact at 1 t_c, independent of dt; worse at the IDORT opacity (Sect. 3.4).
  - **photon bubble, dusty wind, Sincell shock:** not feasible (Sect. 3.6).
- **New defect found.** A hot cell that is optically thick per cell does not emit into a neighbouring thin cell. In the formal solution J_sc ∝ 1/τ_cell, and the moment E reaches only 0.38 of the exact value (Sect. 4). This matters for unresolved photospheres, the IDORT bubbles and wall re-emission in the crooked pipe.

**Spherical counterparts.** The spherical vet_gd counterparts of these tests (homogeneous sphere in 1-D spherical polar, spherical diffusion, etc.) are run separately by the **sptests agent** in /raven/ptmp/jinma/sptests_1011. This note covers only the Cartesian vet_sc path.

## 0. Code, binaries, arms

**Binaries.**
- **Code.** Branch `prodcand-1011`, worktree /raven/u/jinma/ATHENAK/wt-prodcand (read-only for me), at **4c033dcf** ("Merge moveq-1011 ... into prodcand-1011"). It contains cbfix (`implicit_hr_face=sf` default), moveq (`implicit_f0_init`, `implicit_rad_inertia`, `implicit_nofb_work`), icrad/G5, deepatm/gdface (`vet_bc_diff=face` default) and thinfix/vetacc.
- **Binary** bin/athena_cpu_4c033dcf, md5 `1f79726459c6e3484c680b578bdeae07` (copy of /raven/ptmp/jinma/prodcand_1011/bin). It was used for every run unless noted.
- **Later binary.** athena_cpu_e66ee669 (md5 39ab667b48f4e15d29bf13e39f8179ce) adds only opt-in sptests options and docs, so per the coordinator the 4c033dcf results stand. It was used only for the opaque-slab scan (Sect. 4).
- **Test-only build** bin/athena_cpu_4c033dcf_cyl24 (md5 d62869b03950dcab88964e66eb56d4fc):
  - built from `git archive 4c033dcf` in src_bub/;
  - **one change**: the `m1_test = cyl` pgen object limit raised from 8 to 24 (rad_m1_tests2.cpp: `m1_cy_p[24][9]`, `M1CylP v[24][9]`, the limit check and the copy loop);
  - needed for the 20 IDORT bubble sources; solver code identical.

**Arms** (arms.sh; the cen/hr arms were dropped by the user's scope change; earlier cen/hr numbers are quoted from the NOTEs):
- **REC** = `closure=vet_sc vet_tensor=full implicit_flux_faces=all implicit_flux_beam=halfrange implicit_flux=blend implicit_blend=idort_f implicit_blend_tau0=1 alpha=1 flo=0.3 fhi=0.6` + `implicit_blend_xthin=0 implicit_thin_switch=fse implicit_thin_corr=hcs implicit_thin_recon=jse implicit_thin_bq=lin` + `time_scheme=hesdirk2 time2_vet_sc=predict time2_vstage=true implicit_vimp=<true where gas_feedback is on and the gas is dynamic, else false>`.
- **REC-be** = the same keys with `time_scheme=be`.
- Defaults not set: `vet_bc_diff=face`, `implicit_hr_face=sf`, `implicit_f0_init=auto`, `implicit_rad_inertia=true`, `implicit_nofb_work=true`.
- **REC-lag** = REC with `time2_vet_sc=lag`. It is run only where predict fails, as the passing reference the coordinator asked for.
- **vimp.** Marshak and the crooked pipe need `implicit_vimp=false` with REC:
  - Marshak's input does not use bicgstab, and vimp is FATAL there; the production default would resolve it off.
  - The pipe gas is static (kinematic, no gas force).

## 1. Summary table

Pass criteria are those of NOTE-littests Sect. 5 and of the per-branch NOTEs; where none existed, mine are marked (mine). "earlier" = the per-branch result quoted from the NOTEs (cen/hr from littests/cbfix; tJe = the REC keys, mostly BE).

| test | source | pass criterion | REC | REC-be | cen (earlier) | hr (earlier) | status |
|---|---|---|---|---|---|---|---|
| crossing beams, ring 4, 64x256, c dt/dx = 256 | IDORT 4.1.1 Fig 2; J21 4.2 | row energy ≥ 0.9 at top, peak ≥ 0.3 E_b at y = 0.65, floor ≤ 1e-3 | **FATAL** at cycle 2 (predict). REC-lag: 1.006 / 0.417 / 9e-19 | 1.006 / 0.417 / 4e-19 (FWHM 2.6 cells) | 0.74 / 0.074 / 1.9e-2 (FAIL) | 1.006 / 0.417 / 8e-11 (cbfix) | REC-be PASS; **REC FAIL (predict)**; REC-lag PASS |
| same, c dt/dx = 0.5 / 32, ring 12, 256x1024 | " | " | 1.006 / 0.417 at every setting (REC at c dt/dx = 0.5 and 32, ring 12 and 256x1024 all complete) | 1.006 / 0.417 | 0.61-0.74 | 1.006 / 0.417 | PASS |
| crossing beams, Gauss 8x32 (production angles) | " | judge against J_sc (littests) | row 0.91 → 0.83, peak 0.072 (J_sc 0.084) | identical | 2.7-5.2 (creates energy) | 0.92 → 0.85, 0.055 | = cbfix tJe; documented SC-shape limit |
| static shadow, ring 12x2 (36x2), 512x256, γ = 1.000001 | IDORT 4.1.4 Fig 5; J21 4.3 | L1 ≤ 2 %, umbra ≤ 1 % of lit, edges ≤ 1 cell | L1 0.011, umbra 0.22 (0.14 %), wedge 165.3, edge 16 cells (= J_sc) | identical | L1 0.082-0.12, umbra 37 (FAIL) | L1 0.012, umbra 1.4 | PASS |
| static shadow, Gauss 8x32 | " | robustness, vs J_sc | L1 0.030, umbra 0.21 | identical | 0.099 | 0.029 | = cbfix |
| thermal equilibrium 1 (T,E) = (1,100) | IDORT 4.3.1; J21 4.1; JSD12 §4 | E, T → 96.795, 3.13663 | 96.79505, 3.13663 | 96.79505, 3.13663 | same | same | PASS |
| thermal equilibrium 2 (100,1) | " | E, T → 145.788, 3.4747; v = 0 | **E 144.347, T 3.4658, v_x = 1.706**. REC-lag: 145.7879, v 2e-7 | 145.78791, 3.47472, v 0 | same as exact | v 0.021 spurious | REC-be PASS; **REC FAIL (predict)**; REC-lag PASS |
| moving equilibrium, v = 3, c = 10 | IDORT 4.3.1 case 3; J21 | code target E 1.1326, v 2.9554, P_tot 3 (NOTE-moveq) | 1.13329 / 2.95512 / 2.99986 | 1.13333 / 2.95510 / 2.99984 | 1.1333 / 2.9551 (moveq) | 1.1261 / 2.9585, leaks 1e-3 | PASS (old code: FAIL, v 3.001) |
| non-LTE atmosphere, ε = 1e-1 ... 1e-8 | IDORT 4.1.3; J21 4.5 | E vs Eddington analytic: L1 ≤ 2 %, E_top within 1 % | L1 1.1 % (ε = 0.1) ... 0.01 %; E_top/ex 0.9994-1.0000; max dev 10 % at τ ≈ 0.3 (model) | L1 1.08 % ... 0.01 %; E_top 0.9978-1.0000 | surface +6-15 % (marginal) | = tJe | PASS |
| static diffusion, ρκ = 4e4, c = 10 | J21 4.4 | analytic Gaussian, L1 ≤ 1 % | L1 0.02 / 0.06 / 0.28 %, peak 1.0001 | 0.10 / 0.11 / 0.32 %, peak 1.0006 | same as BE | same | PASS |
| dynamic diffusion, gas_feedback = true (ddf), t = 16 | IDORT 4.2.1; J21 4.4 | centroid at v t; peak within 5 % (mine) | peak 0.970, xc - vt = -1e-4, D_eff/D 0.984, L1 0.024 | peak 0.875, xc +1e-4, D_eff/D 1.233, L1 0.115 | 0.875 (BE) | 0.875 (BE) | REC PASS; REC-be position PASS, peak -13 % (BE time error, as before) |
| dynamic diffusion, gas_feedback = false (ddn) | moveq | centroid at v t (old: 4/3 v) | peak 0.998, xc 0, D_eff/D 1.001 | peak 0.930, xc 0, D_eff/D 1.16 | old: +5.33 drift | - | PASS (position); REC-be peak -7 % |
| linear waves, c = 10, P ∈ {0.01, 10}, τ ∈ {0.01, 1, 100}, N 16-256 | IDORT 4.3.2 | L1(ρ) ∝ N^-2; ω matches | 2nd order in all 6 cases (last orders 1.99, 2.03, 2.12, 2.16, 1.43, 2.02); P = 10, τ = 1 L1 2.6e-4 at 256 (tJe earlier 3.8e-4) | 2nd order at P = 0.01; **1st order at P = 10** (BE in time, dt ∝ dx) | 2nd order | FAIL at P = 10 (xthin) | REC PASS; REC-be 1st order = expected BE time error |
| Marshak (T6 gate, S_N reference, c dt/dx = 1 / 10) | JSD12 5.1 | L1(E) ≤ 0.02 | 0.0139 / 0.0141 | 0.0138 / 0.0140 | 0.0135 / 0.0136 | 0.0138 / 0.0140 | PASS |
| irradiated cloud (RHD), ring 12, dt 3e-4, 533 steps | IDORT 4.3.4 Fig 14; J21 Fig 5 | qualitative | x_c(ρ>3) 0.0088 / 0.0157 / 0.0276 at t = 0.06 / 0.10 / 0.16 | identical | 0.003 / 0.005 / 0.009 (fills umbra) | 0.009 / 0.016 / 0.027 | reproduced (= earlier tJe) |
| fixed G5 (2-D pulse, κ 12.8 ... 12800) | g5fix gate | d σ²/dt / 2D within ±2 % | 0.99957 / 0.99987 / 1.00000 / 1.00000 (vimp off) | 0.99958 / 0.99997 / 1.00000 / 1.00000 | 0.99961-1.0 | 0.99975-1.00017 | PASS (= g5fix to 5 digits) |
| **radiation shocks M = 1.2, 2, 3** (new) | IDORT 4.3.3; J21 4.8; JSD12 6.2 | L1 ≤ 2 % in ρ, T_gas, T_rad vs the IDORT transport reference after 3 flow crossings (mine) | 512: L1 ≤ 0.0012 / 0.0011 / 0.0014; 1024 (M2, M3): ≤ 0.0017; precursor -2 ... -9 %, spike ≤ +0.9 % | 512: ≤ 0.0014 | - | - | **PASS** |
| **homogeneous sphere → cylinder** (new, Cartesian analogue) | IDORT 4.1.2; J21 4.6 | E, rF outside within 10 %, f = 1/3 inside (mine) | 256: E L1 out 12 %, F 4.4 %, f_in 0.3344; 512: E 6.1 %, F 5.5 % | 256: 12 % / 4.2 % | - | - | PASS at 512^2, marginal at 256^2 (limit = J_sc, 27 % / 17 % low outside) |
| **crooked pipe** (new) | JSD12 5.2 (Gentile 2001) | no Cartesian reference; dt self-consistency (mine): probes ≤ 5 % between c dt/dx 22 and 220, large dt converges at late times | probes at t = 9.4: 0.468 / 0.427 / 0.385 / 0.343 / 0.276 keV; c22 vs c220 ≤ 1.2 % at t = 1; c2200 vs c220 ≤ 0.8 % for t ≥ 3 | ≤ 0.5 % from REC for t ≥ 1 | - | - | reproduced (qualitative); wall re-emission affected by Sect. 4 |
| **multi-source bubbles** (new) | IDORT 4.2.2 Fig 7 | front at R0 + ct; L1(E) vs exact free streaming ≤ 20 % (mine) | L1 0.43/0.67/0.80 at t = 0.05/0.25/1 t_c (IDORT ρκ); 0.35/0.59/0.63 at ρκ 1e3; lit area at 0.25 t_c 0.23 vs 0.93 | same as REC (1e3) | - | - | **FAIL** (fronts lag the light cone, energy 0.2-0.5 of exact; not dt-dependent) |
| photon bubble / dusty wind / Sincell shock | JSD12 6.4; J21 4.9; JSD12 6.3 | - | - | - | - | - | not feasible (3.6) |

**Regressions against the per-branch results: none.**
- Every REC-be number equals the per-branch tJe number, up to the moveq inertia terms (Marshak 0.0139 → 0.0138) and the face boundary.
- REC improves the BE time error wherever it applies: ddf peak 0.97 instead of 0.875, linear waves P = 10.
- The two REC failures are not regressions, because the same behaviour is on every per-branch binary.
- Bookkeeping: the run.log files of dd*, cloud and two Marshak runs were truncated by an edit to run1.sh while they ran. The athena runs themselves completed (time.txt, all dumps to tlim), and their rc.txt says so.

## 2. Reproduced tests: details

**Crossing beams** (fig/suite_cb.png). REC-be, REC-lag, REC at c dt/dx = 0.5 and 32, ring 12, and 256x1024 all give the cbfix result to the digit:
- row energy 1.006 at every height (1.128 in the top row next to x1max);
- peak 0.417 at y = 0.65, FWHM 2.6 cells (the 2-cell zigzag of NOTE-cbfix);
- floor ≤ 1e-16.

**REC at c dt/dx = 256 (J21's dt) is FATAL at cycle 2.**
- hesdirk2 stage 1 is "NOT ADMISSIBLE" (min E < 0). The backward-Euler redo then diverges: Picard resid 2e25, `implicit_resid_fatal`.
- With `time2_vet_sc=lag` the run completes (Picard 3.1). It also completes without vstage and lag.
- In the empty box hesdirk2 stage 1 is inadmissible on most steps: 16/20 at c dt/dx = 256, 1099/1100 at 0.5, 23/24 at 32. So in vacuum transport REC runs as BE in practice.

**Static shadow** (fig/suite_shd_r12g2.png).
- REC and REC-be are identical to cbfix/littests tJe: L1 0.011, lit 158.2, wing 79.3, wedge 165.3 (exact 157.6), umbra 0.22 at x = 0.15, edge 16 cells (= J_sc).
- Gauss 8x32: L1 0.030.

**Thermal equilibrium** (fig/suite_thermal_eq.png, fig/suite_eq2_predict.png).
- Case 1 is exact in both arms.
- **Case 2 under REC:** the mean F_1 reaches 4.4e3 = 0.30 cE at t = 0.015 in a uniform periodic box.
  - The gas is accelerated to v_x = 1.71 (about Mach 0.7), and E ends at 144.35 instead of 145.79.
  - Turning off vimp, vstage or `implicit_rad_inertia` changes nothing.
  - `time2_vet_sc=lag` gives the exact state, |v| 2e-7.
  - Same numbers on 42aa6bbc, 9e1f2686 and mq3.

**Moving equilibrium:** REC and REC-be reproduce NOTE-moveq (E 1.1333, v 2.9551, P_tot 2.99985). The old code failed (v 3.001).

**non-LTE / static diffusion / Marshak / G5:** as in the table. These are local or diffusive tests; the REC/REC-be differences are at the 1e-3 level (fig/suite_nlte.png, fig/suite_diffusion.png).

**Dynamic diffusion** (gas_feedback on): REC fixes the BE peak deficit (0.970 against 0.875). The NOTE-moveq h2 value is 0.970.

**Linear waves** (ana/rw_table_suite.txt, fig/suite_rw.png).
- REC is 2nd order in every case and matches ω_ref to < 1e-3. At P = 10, τ = 1 it is better than littests tJe (2.6e-4 against 3.8e-4 at N = 256).
- REC-be is 1st order at P = 10. That is BE time error with nt = 2N, not a spatial defect: at P = 0.01 the two arms agree to 3 digits.

**Irradiated cloud** (fig/suite_cloud.png): REC = REC-be = littests tJe morphology and centroids.

## 3. Formerly "not done" tests

### 3.1 Radiation shocks (IDORT 4.3.3, J21 4.8). PASS. Figure: fig/suite_shock.png

**Set-up.**
- Lowrie-Edwards parameters (P0 = 1e-4, σ_a = 1e6, κ = 1, i.e. C = 1732, σ_t = 577), M = 1.2, 2, 3. This is exactly the IDORT set (a = 1e-4, c = 1.73e3, σ = 577.4).
- `radshock` pgen in the cgs scaling of tests_m1/t7_radshock.py, which also writes the IC (the LE08 non-equilibrium-diffusion profile; inp/rs_m*.athinput).
- 512 cells (1024 for REC at M = 2, 3) on 0.0152-0.0157 cm, nx2 = 4.
- **Boundary.** The "implicit Dirichlet/efix" set-up asked for is replaced by **Marshak at both ends with bath = a T^4 of each equilibrium**.
  - This is the runs_3z_validate finding: efix breaks the outflow end.
  - The hydro user BC holds the far-field states.
- **Run length:** 3 flow-crossing times L/v0 (1.3e-9 s at M = 2), at the hydro CFL 0.4 (c dt/dx about 240). That is 5800-10800 steps (about 23000 at 1024), 5-15 min per run.

**Reference.**
- The **IDORT/J21 transport solutions** from the IDORT Zenodo package (record 15032138, IDORT_test_mini.zip, md5 8505879774d008697fb44a23b4e3d305, `shock_1d/Rad_shock_sol/M{1p2,2,3}`; 4096-cell S_N steady profiles, the "Ferguson+17" reference).
- Profiles are aligned on the density jump and compared over the code domain (ana/shockcmp.py).

| run | L1 ρ / T_gas / T_rad at t = 3 L/v0 | precursor (1 % level) vs ref | T_gas spike vs ref | IC (LE08 diffusion) L1 |
|---|---|---|---|---|
| M1.2 REC / REC-be, 512 | 4e-4 / 4e-4 / 3e-4 (both) | -9.3 % | +0.02 % | 1.5e-3 / 1.0e-3 / 1.9e-3 |
| M2 REC, 512 / 1024 | 7e-4 / 1.1e-3 / 1.4e-3; 4e-4 / 6e-4 / 7e-4 | -2.3 % / -8.9 % | +0.8 % / +0.6 % | 4.1e-3 / 6.1e-3 / 1.1e-2 |
| M2 REC-be, 512 | 7e-4 / 9e-4 / 1.3e-3 | -2.3 % | +0.7 % | |
| M3 REC, 512 / 1024 | 1.2e-3 / 9e-4 / 8e-4; 1.0e-3 / 1.7e-3 / 1.7e-3 | -2.9 % / -0.9 % | +0.9 % / +0.2 % | 9.9e-3 / 2.4e-2 / 2.6e-2 |
| M3 REC-be, 512 | 1.0e-3 / 1.4e-3 / 1.1e-3 | -2.9 % | +0.4 % | |

**Results.**
- The VET solution relaxes away from the diffusion IC (precursor -12 ... -25 %) onto the transport reference, then holds steady with a shift ≤ 1.5 cells.
- Against the LE08 diffusion reference (t7_radshock.py), M = 3 reads FAIL (T_rad L1 2.1 %). The cause is the precursor, which is 15 % longer than in diffusion. That is the expected transport-vs-diffusion difference, and the transport reference is the right one to judge against.
- Solver: Picard 2.0-3.7, 0 NON-CONVERGED, 0 hesdirk2 fallbacks.

### 3.2 Homogeneous sphere → Cartesian homogeneous cylinder (IDORT 4.1.2 / J21 4.6). Figure: fig/suite_hcyl.png

**Why a cylinder.** vet_sc is Cartesian, and the existing `cyl` pgen is 2-D, so a sphere is not available here. The 1-D spherical test belongs to the sptests agent.

**Set-up** (inp/hcyl.athinput):
- Infinite cylinder R0 = 1, ρκ_a = 10, T = 1 held (γ = 1.000001, gas_feedback off), vacuum outside (ρ 1e-10).
- Box [-4,4]^2 at 256^2 (τ_cell = 0.31) and 512^2, x1 vacuum Marshak, x2 periodic.
- c = 100, c dt/dx = 100, 64 steps (128 at 512^2), to steady state.

**Exact reference** (ana/hcyl.py):
- I = B(1 - e^{-τ}) summed over the periodic images, by 3-D ray quadrature (64 μ x 4096 φ).
- The analytic check passes: E = 1 and f = 1/3 at the centre, E = 0.5 and F = E/2 at the surface.

| run | E inside L1 | E(R0)/ex | E outside L1 | F outside L1 | r F (2-3.5) / ex | f_xx at r = 0.5 (ex 0.3331) | f_xx at r = 3 (ex 0.297) | J_sc outside L1 |
|---|---|---|---|---|---|---|---|---|
| REC 256^2 | 9e-4 | 0.967 | 0.118 | 0.044 | 1.034 | 0.3344 | 0.367 | 0.266 |
| REC-be 256^2 | 8e-4 | 0.965 | 0.120 | 0.042 | 1.032 | 0.3344 | 0.367 | - |
| REC, nphi 128 | 9e-4 | 0.967 | 0.117 | 0.043 | 1.035 | 0.3344 | 0.366 | 0.265 |
| REC 512^2 | 1.4e-3 | 0.997 | 0.061 | 0.055 | 1.043 | 0.3331 | 0.339 | 0.171 |

**Reading.**
- Inside, the diffusion limit is exact (f = 1/3, E = 1).
- Outside, the moment E is 8-15 % low at 256^2. The limiting error is the formal solution: J_sc is 27 % low, because bilinear short characteristics smear the emerging field.
- The error converges at first order: 12 % → 6 % in E, J_sc 27 % → 17 %.
- Not an angle effect: nphi 128 gives the same.

### 3.3 Crooked pipe (JSD12 5.2, Gentile 2001 geometry in Cartesian). Reproduced. Figure: fig/suite_pipe.png

**Set-up** (inp/pipe.athinput):
- (0,7)x(-2,2) at 512x256. The 7 opaque walls of JSD12 (ρ 10, σ_a 2000, τ_cell 27) are `cyl` rectangles; the pipe is ρ 0.01, σ_a 0.2 (powerlaw κ = 200 (ρ/10)^{1/3}).
- c = 300. Units are cm, 1e-8 s and 0.5 keV; c_v = 1e15 erg/g/keV, giving γ - 1 = 1/58.32. T0 = 0.05 keV in equilibrium.
- **Source.** The isotropic 0.5 keV source is an `implicit_ebath_x1min = 1` bath on the whole x1min face. The walls cover x1min for |y| > 0.5, so this is JSD12's |y| < 0.5 source; the wall faces only warm their first cell.
- **The patch source as a hot held cell fails:** see Sect. 4. A ρ = 1e6 (τ_cell 1e8) held source emitted nothing into the pipe.
- x2 periodic: the top and bottom rows are wall.
- Static gas: kinematic + advect + `dbg_gas_force(_trans)=false`, max|v| = 0 verified. This is JSD12's "velocities kept zero".
- Probes at Gentile's (0.25,0), (2.75,0), (3.5,1.25), (4.25,0), (6.75,0).

**Results.**
- REC at c dt/dx = 220 (dt = 1e-2), T at the probes (keV):

  | t | probe 1 | 2 | 3 | 4 | 5 |
  |---|---|---|---|---|---|
  | 0.8 | 0.450 | 0.385 | 0.322 | 0.267 | 0.182 |
  | 9.4 | 0.468 | 0.427 | 0.385 | 0.343 | 0.276 |

  - The front passes the probes in order: T ≥ 0.25 keV at t = 0.1 / 0.1 / 0.3 / 0.7 / 4.0.
  - Walls stay at 0.05 keV.
- **dt consistency.**
  - c dt/dx 22 vs 220: max probe difference 19 % at t = 0.3 (front arrival), 5 % at 0.5, 1.2 % at 1.
  - c dt/dx 2200 (dt = 0.1) vs 220: misses the early front (factor 4 at t = 0.3), then 0.8 % at t = 3 and 0.4 % at t = 10. This is JSD12's own large-dt behaviour (their dashed line).
  - REC-be vs REC: ≤ 2.7 % at t = 0.5 and ≤ 0.1 % from t = 3.
- **No quantitative reference exists in Cartesian geometry.** JSD12 say themselves that they cannot reproduce Gentile's RZ numbers, so the status is "reproduced qualitatively".
- Wall re-emission into the pipe is suppressed by the τ_cell = 27 wall surface (Sect. 4). With a resolved wall surface the late-time probe temperatures would be somewhat higher.
- Cost: 685-2500 s per run on 18 threads; Picard 3.0-3.6.
- The Caltech TASK (TASK-2026-10-11-caltech-crooked-pipe.md) was withdrawn because Raven finished first.

### 3.4 Multi-source bubbles (IDORT 4.2.2). Figure: fig/suite_bub.png

**Set-up** (inp/bub.athinput, binary 4c033dcf_cyl24):
- The 20 sources are the IDORT ones: create.py's logistic-map recipe, positions in [-0.4, 0.4]^2, r ≤ 0.03, T ≤ 10.
- They are opaque discs held at T (`m1_test = cyl`, user_srcs; gas static: kinematic + advect, γ = 1.000001, gas_feedback off). The vacuum is ρ 1e-12, κ_a ρ ≈ 0.
- Periodic unit box (the IDORT param.txt has PeriodicBoundariesOn = 1), 512^2. c = 10, dt = 1e-3 (c dt/dx = 5.12, IDORT's dt_max), 100 steps to t = 0.1 = 1 t_c.
- Angles: ring 100 (in-plane, like IDORT's 200 2-D rays) and Gauss 8x32.
- Opacities: disc ρκ = 1e5 (IDORT; τ_cell 195) and 1e3 (τ_cell 1.95).
- Exact reference (ana/bub.py): free streaming from the discs' surfaces, I = T_i^4 of the first disc hit when it lies inside the light cone. "E2" uses in-plane rays (for the ring) and "E3" uses 3-D directions (for the Gauss set). It is checked against the code's own steady J_sc, which is 0.94-0.96 of E2(t = 0.1) at 128^2.

**Metrics** at t = 0.005 / 0.025 / 0.1: L1(E)/ΣE_ex; total E/ΣE_ex; the "lit" fraction of the box with E > 1 (code / exact).

| run | L1 | total E / exact | lit fraction (code / exact) | wall |
|---|---|---|---|---|
| REC ring100, ρκ 1e5 (IDORT) | 0.43 / 0.67 / 0.80 | 0.57 / 0.33 / 0.21 | 0.10/0.17, 0.23/0.93, 0.58/1.00 | 236 s |
| REC ring100, ρκ 1e3 | 0.35 / 0.59 / 0.63 | 0.66 / 0.46 / 0.42 | 0.10/0.17, 0.24/0.93, 0.60/1.00 | 188 s |
| REC-be ring100, ρκ 1e3 | 0.35 / 0.59 / 0.63 | 0.66 / 0.46 / 0.42 | same | 136 s |
| REC Gauss 8x32, ρκ 1e3 (vs E3) | 0.32 / 0.55 / 0.64 | 0.69 / 0.48 / 0.42 | 0.09/0.16, 0.21/0.93, 0.53/1.00 | 641 s |
| 128^2, ring 100, ρκ 250 (τ_cell 1.95): REC at c dt/dx = 0.32 / 1.28 / 5.12; REC-be 0.32 | 0.33-0.29 / 0.55-0.51 / 0.59 | 0.69-0.80 / 0.49-0.52 / 0.43 | 0.15-0.28/0.17, 0.39-0.52/0.94, 0.92-0.94/1.00 | - |

**Verdict: FAIL** against the IDORT figure. The bubbles do not fill the box at R0 + ct.
- **Front.** Along a ray from an isolated disc at 128^2 the E front sits at about 0.75-0.8 ct, smeared over many cells. The moment E stays at 0.45 of J_sc in total at t = 1 t_c.
- **It is not the time step.** c dt/dx 0.32 → 5.12 changes nothing at t ≥ 0.025, and REC = REC-be. It is not the angle set either (ring 100 ≈ Gauss).
- **At the IDORT opacity** (τ_cell 195) the deficit is worse: 0.21 against 0.42. That is the opaque-emitter defect of Sect. 4, which suppresses the discs' luminosity.
- **Likely cause** (not verified): the VET and the half-range ratios come from a **steady** formal solution. Ahead of a front they describe the steady field of all discs, not the radiation that has arrived. In IDORT each ray is time-dependent. This is the "moment methods fail" point the IDORT test was designed to make.
- The appearance is otherwise sensible: bubbles merge, there is no spurious shadowing, and the ring-100 run shows DOM-like ray streaks at t = 1 t_c at ρκ = 1e5.

### 3.5 Inputs from the IDORT package

- The IDORT Zenodo record 15032138 is linked from the arXiv comments.
- `IDORT_test_mini.zip` (36.5 MB, md5 checked) was downloaded to ext/idort_mini and treated as data. Only these were used:
  - `bubble_2d/create.py` (the deterministic logistic-map recipe for the 20 sources: positions, radii, temperatures; reproduced in ana/bub.py → ana/bub_sources.json);
  - `shock_1d/Rad_shock_sol/*` (shock references).
- The JSD12 LaTeX (arXiv e-print) was fetched for the crooked-pipe parameters (ext/jsd12).

### 3.6 Not feasible here

- **Photon bubble (JSD12 6.4).**
  - It needs MHD with a radiation-supported stratified atmosphere and the Blaes-Socrates linear growth rate.
  - There is no such pgen. The `m1_test` pgens that accept `<mhd>` are only the uniform/beam/pulse/wave ones, and the implicit vet_sc + MHD coupling has not been validated.
  - Writing the pgen plus an analytic dispersion check is a project of its own, not a test run.
- **Dusty wind (J21 4.9).**
  - It needs a new pgen (gravity, dust opacity, wind launching), at 512x1024 for 2e4 steps.
  - Cost estimate from the bubble timing (about 28 s per 512^2 step with 100 rays on 16 threads): about 1e6 s of CPU for one run. That is out of scope without a GPU vet_sc path for this geometry.
- **Sincell VET shock (JSD12 6.3):** no pgen and no reference profile in hand. The LE/IDORT shocks of 3.1 cover the same physics.
- **JSD12 ±14° shadow:** not representable on a ring that contains ±14° without a ray along x2 (NOTE-littests). The ±15° IDORT shadow is the one run.

## 4. New defect: a per-cell-opaque hot cell does not emit into a thin neighbour

**Where it was found.** It turned up while building the crooked-pipe source, and it is confirmed by a clean 1-D slab scan.
- Inputs: inp/slab.athinput, runs/slab, binary e66ee669, REC.
- Set-up: a held T = 1 slab of width 0.1 (12.8 cells) next to vacuum. At x = 0.5 the exact field is E = 0.5 aT^4 and F = c aT^4/4.

| τ_cell of the slab | E/E_ex | F/F_ex | J_sc/E_ex | H_sc/F_ex |
|---|---|---|---|---|
| 0.1 | 0.996 | 1.058 | 0.894 | 0.950 |
| 1 | 0.789 | 0.953 | 0.557 | 0.672 |
| 10 | 0.441 | 0.588 | 0.100 | 0.133 |
| 100 | 0.381 | 0.509 | 0.010 | 0.013 |
| 1e4 | 0.375 | 0.500 | 1e-4 | 1e-4 |

**Cause** (rad_m1_vet.cpp, the SC segment, around l. 930).
- The segment from the upwind point to the cell centre uses the trapezoidal dτ = (χ_u + χ_0)dx/(2μ) and a source that is linear in τ between S_u and S_0.
- When the opacity jumps, all of dτ sits in the opaque half, but the linear source puts S_0 (vacuum, 0) at the downwind end.
- The emergent intensity is therefore S_0 + (S_u - S_0)(g - e^{-dτ}) ≈ S_u/dτ instead of S_u.
- The moment solve then keeps only about half of the flux (E 0.38, F 0.50 of exact for τ_cell ≥ 100).

**Consequences.**
- Any unresolved hot photosphere next to thin gas is affected.
- In the pipe test, a ρ = 1e6 source cell (τ_cell 1e8) emitted nothing: pipe E stayed at 1e-4 and J_sc at 1e-4.
- In the IDORT bubbles at the paper's ρκ = 1e5 (τ_cell 195 at 512^2), see 3.4.

**Suggested fix (not done).** Integrate the segment as two half-cells, each with its own χ and S:
`I = I_u e^{-(τ_a+τ_b)} + S_u (1 - e^{-τ_a}) e^{-τ_b} + S_0 (1 - e^{-τ_b})`, with τ_a = χ_u dx/(2μ) and τ_b = χ_0 dx/(2μ).
Alternatively, switch to this form only where χ_u/χ_0 is large.

## 5. Files

- **Scripts.**
  - run1.sh, run1b.sh (bubble/test binary, BIN overridable), arms.sh, pool.sbatch / poolb.sbatch (general, 72 cores), rw.sbatch, g5/ (run_g5.sh, g5.sbatch), build_bub.sh.
  - Job lists: jobs_s1, jobs_sh, jobs_hc, jobs_pipe, jobs_bub.
- **Inputs** (inp/):
  - littests inputs plus eq3t, ddn, ddf from moveq;
  - new: rs_m{1.2,2,3}, hcyl, pipe, bub, slab.
- **Runs** (runs/<test>/<arm>_<tag>): runs/rw (linear waves), g5/run (G5).
  - runs_aborted1/ holds superseded attempts: the cen/hr runs cancelled by the scope change, smoke runs, the wrong-binary bubble attempt, the hydro-dt-limited bubble attempt and the failed-vimp Marshak runs.
- **Analysis** (ana/):
  - littests/cbfix/moveq scripts: cbs.py, shdsum.py, s2.py, mq/eqts.py, mq/ddsum.sh, rwtab.py;
  - new: shockcmp.py, hcyl.py, pipe.py, bub.py, figs.py, fig_cb.py, fig_shd.py, fig_1d.py, fig_cloud.py;
  - exact caches: hcyl_exact_{256,512}.npz, bub_exact_512_8.npz;
  - tables: rw_table_suite.txt.
- **Figures** (fig/, white background):
  - suite_cb.png, suite_shd_r12g2.png, suite_thermal_eq.png, suite_eq2_predict.png, suite_nlte.png, suite_diffusion.png, suite_rw.png, suite_cloud.png;
  - suite_shock.png, suite_hcyl.png, suite_pipe.png, suite_bub.png.
- **Handover for the pipe:** wt-bsg docs/handover/suite-pipe-1011/ (committed locally as 91efa580), plus the TASK file (withdrawn).
