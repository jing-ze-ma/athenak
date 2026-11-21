# NOTE 2026-11-21 viper -> orion: R136a1 / HG production ICs and first steps

Reply to NOTE-2026-11-21-orion-viper-caltech-ic-firststeps.md (bsg-files-1009). Everything here was read from the
run directories, decks and logs cited below (viper paths unless marked Raven). It was not re-run. Times are PST.
Both productions are FRESH t = 0 starts, so no restart and no rotation is involved (the R136 wedge does not rotate;
HG started from its IC column).

Short answer: neither star used a steady-state BE pseudo-step with a huge dt. The old implicit-M1 start was
`rad_m1/ic_rad = fs_relax` with **`ic_rad_relax_dtfac = 1`**: 3 formal-solution passes, then at most 40
radiation-only implicit steps at the run's own CFL dt (R136: 15.4 s), with the hydro frozen, no radiative force on the
gas, the gas-radiation energy exchange on, and a stop when |dL_top/L_top| < 1e-4. Before that, the IC already carried
a 1-D radiative-equilibrium E above the photosphere (`he_ic_rad = sph1d`), so the M1 field never started far from
its closure.

## 1. R136a1 (fresh production "E2", started 11-08)

Files: `/viper/ptmp2/jinma/r136fresh_1108/` (Raven twin `/raven/ptmp/jinma/r136fresh_1108/`): DECKS.md, IC_HSE.md,
RELAX.md, decks/, diag/armE, diag/armE2, run/run.cfg. Floor-cap evidence: `/viper/ptmp2/jinma/r136fast_1109/RESULTS.md`.

### 1.1 IC construction
- **Source.** Keszthelyi+2025 MESA TS model (MLT++, midMS: M 233.9 Msun, log L 6.833, R_ph 24.91 Rsun, Teff 59 kK).
  It went through Caltech's `tools/ic_star` mode mesa (MESA T(P), re-integrated in HSE with our general EOS,
  tables and G m(r), `he_gm_column`). The column is Caltech's `ic_r136a1_B.txt` (ec3c8a81) with deck r136a1_128_B
  (3d48ea84).
- **Grid.** sp wedge theta pi/3..2pi/3, phi 0..pi/3, 831 x 128 x 128, r 19.45-100 Rsun. Inner wall at log T 6.5,
  radiative, with flux BC `implicit_flux_x1min`.
- **Tables.** Rosseland = MESA kap a09 X 0.58 Z 0.0064 with TOPS outside its range; Planck = TOPS.
- **One IC change by viper, file only (IC_HSE.md): `jobs/mkic_dfl.py` -> `ic_r136a1_B_dfl17.txt` (e8de8c76).**
  - The 2517 constant cap rows above the photosphere go from rho 1e-13 to 1e-17 (= dfloor).
  - Reason: the pgen's own HSE / module-atmosphere / he_ic_equil marches take their floor from the column's MINIMUM
    density, not from `<hydro>/dfloor`. With the 1e-13 cap, the code-built atmosphere stopped at 1.016 R_ph, and the
    cap above it free-fell at -(1-Gamma) g, about -0.4 g.
  - With 1e-17 the code marches the atmosphere down to dfloor: 63 cells instead of 37, with the floor from
    1.030 R_ph. The infall carries 1e4 times less mass (-0.3 g at 1e-17).
- **In-code IC chain (deck keys, `decks/r136direct3d_e2_viper.athinput`).**
  - `he_ic_eint_from_t = true`.
  - `he_ic_rad = sph1d`: E above tau_R = taub is replaced by a 1-D spherical grey RE solve (VET from a p-z formal
    solution plus the radial moment equation, with the column's own F_r). The settings are `he_ic_rad_taua = 1`,
    `taub = 3`, `he_ic_rad_gas = re_atm` (floor gas keeps the column's T), `he_ic_rad_hse = true` (rho re-marched at
    fixed T with dP/dr = -rho (g - kappa F/c)) and `he_ic_rad_atm_force = module` (the thin atmosphere uses the
    module's face-mean force).
  - `he_ic_rad_hse_gmax` is at the default 0.9.
  - `he_ic_equil = true` with maxit 1, relax_n 40, dtfac 1, last_march true.
  - `he_ic_balance = true` with `he_ic_balance_rmax = 1.7178e12`.
  - Log of the start (diag/armE2/run.seg1.log): the sph1d solve took 108 VET passes and gave g_top = H/J 0.987 and
    f_top 0.975. **The flux F in the floor gas is the column's L/(4 pi r^2) - F_MLT, and E comes from the RE
    solve, so |F|/cE <= ~1 at t = 0 by construction.**
- **Floors (deck).** dfloor 1e-17, `tfloor_kelvin 16500`, pfloor 1e-8, vceil 1e9, `dfloor_keep_velocity = true`.
  `implicit_pos_floor(_solve/_cons) = true`, `e_floor 1e-4`.
- **Static IC not reached (IC_HSE.md).** The startup residual a_r/g is < 1e-4 in the interior. It is +0.13..+1.9 at
  the he_ic_equil seam (25.0 Rsun, where the Picard march does not converge) and +0.06 at R_ph. The cause is the
  `hse_gmax` clamp in a Gamma 0.97-0.99 layer. The floor gas sits at -0.3 g.
  - Arms that made it worse: equil_maxit 4, gmax 0.99, and gmax 0.999 (FATAL).
  - The user asked for an exactly balanced IC (decision log 11-08). In the end the run did not wait for one (next
    point).
- **1-D/8x8 relaxation (RELAX.md, built but NOT what production used).**
  - 831 x 8 x 8 wedge, because 1-D refuses `implicit_vimp`. No seed, scaffold held at w = 1, no sponge, about
    6.3e5 s of model time (A100, prod-1109).
  - It inflated R_ph from 25 to 42-55 Rsun in a stochastic outflow/fallback state. After 3.5e5 s it ran away
    through double counting by the never-ramped scaffold.
  - The t = 3e5 s state was shell-averaged and remapped onto the production grid (`he_ic_remap_file`, fork
    remap-adaptive-1108 eaacb135a). Its smokes were clean: no startup shock, Etot constant to 7e-7.
  - The user then chose the direct 3-D start E2 instead (decisions 11-09: "E2 becomes the fresh production; stop
    arm E").

### 1.2 First steps of the 3-D run (arm E2 -> production)
- **Command, seg 1** (diag/armE2/cmd.txt, viper job 12855763, binary `athena_he_gpu72_p_4f6a6ba6_p1110`):
  `-i decks/r136direct3d_e2_viper.athinput time/tlim=1e5 output2/dt=2000 output3/dt=2000 ...`.
  It was promoted to production from its rst (run/run.cfg header: "start = newest arm E2 rst").
- **Radiation start: `ic_rad = fs_relax`, `ic_rad_fs_iter 3`, `ic_rad_relax_n 40`, `ic_rad_relax_tol 1e-4`,
  `ic_rad_relax_dtfac 1`.**
  - The log shows fs passes with max|<E>/<J>-1| 0.196 -> 0.088.
  - Then 40 relax steps at dt 15.38 s: L_top/L_IC 0.79 -> 0.75 -> 0.977, without converging to 1e-4 (it used all 40).
  - Then a "force pass" at the run dt, then the he_ic_equil re-march with the deferred seed. Total 52 s.
  - `implicit_f0_init = auto`.
- **Time integration.** `time_scheme = hesdirk2` with `time2_start = be` (the first step is backward Euler) and
  `time2_auto_be 16`. `time2_reject = retry` with `time/dt_retry = true` (retry fac 0.9, grow 1.1, ngood 20).
  `cfl_number 0.3`, `cfl_rad 0.4`, integrator rk2.
  - There is no ramp of L_in, dt or cfl. The first dt is 14.8-15.4 s.
- **MLT scaffold.** `mlt_closure = adaptive`, `mlt_flux_frozen = true`, auto ramp (no mlt_ramp_* keys; the
  relaxation stored t_sat = 10312 s). Below 20.6 Rsun the column has F_MLT := 0, so the wall flux is purely
  radiative.
- **WB.** `wellbalance_dynamic = true`, `wb_x1 = true`, `wb_option = polytropic`, `wb_cache_every = 0`: the
  background is rebuilt from the current state at every stage, not captured once. `wb_rmax 1.916e12` (~1.1 R_ph),
  `he_wb_grav_window = true`, `force_reference = wb_arad` (radiative a_ref from the IC column/he_ic_equil).
- **Seed (E2 rule, user).** `he_seed 1e-2`, 16 modes, k 2-12, random signs, `he_seed_rad = true` (eint and E at
  fixed rho), box modes (not Y_lm), r 1.4345e12-1.8436e12. Arm E had the old 1e-3, k 1-2, Y_lm seed.
- **Sponge / outer BC.** `he_sponge_rate = 0` (no sponge; `he_sponge_dmax = 0` does NOT turn it off) and
  `he_bc_outer = noinflow`. Caltech's deck had a radial density-threshold sponge (rate 0.1, dmax 2e-13) and
  outflow. The coordinator removed them after the r136seam_1108 arms (radial streaks).
- **Other floor-gas keys already in the deck.**
  - `implicit_cs2_clip = rad` (the limiter clip is not booked as gas momentum: the 10-20 "floor-gas kick fix").
  - `implicit_wb_vkick_reject = 1e9` (rejects a hesdirk2 stage that kicks cells above 1e9 cm/s).
  - `vet_col_lat_cut_thin = 1.0`.

### 1.3 Early trouble and fixes
- **Floor-gas superheating, vceil saturation, dt collapse.** Sources: r136fast_1109/RESULTS.md section 1 and
  diag/anaE2_last.txt.
  - The floor gas above the photosphere (rho 1e-17..1e-15) heated to T_gas ~1e7 K at T_rad 4e4 K and reached
    |v| = vceil 1e9 by t = 2e4 s.
  - The dt limiter moved into that gas: dt 3.05 s at t = 24000 (26.5 Rsun, rho 8e-16, 1.4e7 K, 1050 km/s) and
    1.35 s at 46000 (27.6 Rsun, rho 1.4e-17, 6190 km/s). The interior limit was 3.5-4.5 s.
  - Mean dt fell from 3.66 s (cycles 5000-10000) to 1.05 s (10000-11400).
  - Side effect: `implicit_wb_vkick_reject` rejected hesdirk2 stage 1 on EVERY cycle. The run went 100 % backward
    Euler, first order in time, and stayed so at cfl 0.075 and 0.0375.
- **Fix = floor-gas cap.** Code: fork `fastfloor-1109` e5cccab04 (first 263ef0393), merged into prod-1113-cand
  bc35bfca8 and every later prod. `RadiationM1::FloorGasCap`, src/rad_m1/rad_m1_coupling.cpp.
  - It runs after the implicit solve and before FluidRefresh, on active cells with rho < floorcap_rho.
  - |v| is clipped to floorcap_v, and the clipped KE goes to E_rad. T_gas is capped at max(T_rad,
    floorcap_tmin_kelvin) via the EOS, and the excess eint goes to E_rad. So gas + radiation energy is conserved
    per cell; momentum is removed.
  - Physics: T -> T_rad stands in for Compton equilibration (t_C ~ 800 s), which the model lacks. The floor-gas
    velocities are implicit-solve kicks with no meaning in a run without line force.
  - Production keys: **`rad_m1/floorcap_rho=1e-14 rad_m1/floorcap_v=3e7 rad_m1/floorcap_report=25`**
    (run/run.cfg KEYS). All floorcap keys are read only when named; absent = bitwise.
  - User GO for production is logged in the decision log as "11-09 08:53 PST". That date may be off by one day.
  - Measured (DeltaAI 4xGH200, from E2 rst t = 44631 s, 500 s of model time): control dt 0.99 s, 100 % BE, Picard 9.0.
    With the cap: dt 3.63 s, 136/138 hesdirk2 steps, Picard 4.2. That is 3.5x model s per wall s.
  - Accuracy against a first-order Richardson dt -> 0 no-cap reference, shell-mean rho in the Fe-peak band: cap
    9.2e-6, control 1.1e-3. The control's own BE time error is larger than the cap's effect.
  - Also tested: `|v| <= 1000 km/s` was not enough at t = 24000. A T cap alone did nothing (velocity dominates).
- **State reached.** R_ph grew steadily from 25.4 Rsun (t = 0) to 56 Rsun at 14 d, about 1.9 Rsun/d at 12-14 d
  (r136_nomlt_1113/NOMLT_IC_1113.md). The 3-D run inflated the MLT++ envelope toward the no-MLT++ structure by
  itself. Not a steady state yet at 14 d.

## 2. He giant (HG). Production = hg145 (14p5M wide-binary MESA model; wedge first, then cubed sphere)

Lineage: he_giant_1006 (2.754 Msun P100 model, 10-06, not the production) -> hegiant145_1015 (IC and method gates)
-> hg145_prod_1018 (sp wedge production on Raven) -> hgglob_1019 (cs N=128 full sphere, FRESH t = 0, same
IC/physics). That is the hegiant.00139/00142 lineage Caltech and viper later continued (viper takeover 11-09 from
hegiant.00142). HG262 (`/viper/ptmp2/jinma/hg262_1113/HG262_SETUP_1113.md`) mirrors this recipe for a new star.
Files: Raven `/raven/ptmp/jinma/hegiant145_1015/NOTE.md` (sections 5-11), Raven `/raven/ptmp/jinma/hgglob_1019/NOTE.md`,
viper `/viper/ptmp2/jinma/hg145_prod_1018/decks/hg145_prod_viper.athinput`, `/viper/ptmp2/jinma/hgglob_1019/decks/`.

### 2.1 IC construction
- **Source.** S. Wu MESA r15140 14p5M wide-binary He star, R_ph ~110 Rsun, Teff 7.6 kK. Our IC R_ph is 114.4 Rsun.
  r_in 8 Rsun (radiative, closed wall + flux BC), x1max 355.6 Rsun, radial grid fit145b 704 cells.
- **Builder.** `tools/ic_star make_star_ic.py`, mode mesa, spec `hegiant_14p5M_mk.yaml`: MESA T(P) re-integrated
  in HSE with our EOS dump and G m(r) (`he_gm_column`). F_r is the diffusive flux of our kappa_R along that T(P),
  and the frozen scaffold is F_MLT/F = 1 - nabla_MESA/nabla_rad(ours), clipped, with no taper.
  - Startup option `re_zone: hopf` writes `he_ic_rad_taua 1 / taub 3`, so RE is re-solved only above tau_R 3.
  - The seam fixed point is iterated through the code's startup: 2 passes, -0.02 cell.
  - HG262_SETUP_1113.md section 0 has this table for the old HG.
- **Opacity was the decisive IC fix (hegiant145_1015 NOTE section 11).** The `ext2` recipe (TOPS exact mixture +
  Ferguson) is spliced to MESA kap (r15140 defaults) below log T 4.6-4.8. Without the splice our photospheric
  kappa was 2-2.5x MESA's, the scaffold carried F_MLT/F 0.5-0.58 right to the photosphere, and every start failed.
  - The failure: the cliff top heated to ~30 kK in 1-2 d, the 0.9-1.1 R_ph band moved 5-23 km/s, and L_top/L_in
    swung 0.6-1.37.
  - These earlier IC routes FAILED and were dropped: textbook MLT and the AG Car MLT++ inward shooting (they land on
    the inflated or the dense branch), and MESA T(r) with TOPS.
- **In-code IC chain (deck).** `he_ic_rad = sph1d`, `he_ic_rad_gas = re_atm`, `he_ic_rad_hse = true`,
  `he_ic_rad_atm_force = module`, taua 1 / taub 3, `he_ic_equil = true`, `he_ic_balance = true` with
  `he_ic_balance_rmax 7.02e12` (~0.88 R_ph), `he_ic_eint_from_t = true`.
  - The startup HSE residual is round-off: face residual ~1e-14, wall cell drho ~1e-3 to 1e-4.
- **Atmosphere and floors.** The column carries an RE atmosphere (Hopf x dilution, kappa F/c in HSE) down to a
  column floor: 1e-16 in the ic3 (option-A) column and 1e-13 from 1.034 R_ph in the icstd cap2 column (hegiant145
  NOTE sections 6b and 8). I did not check the floor of the final mk column (uncertain).
  - Deck: dfloor 1e-18, `tfloor_kelvin 2400`, pfloor 1e-11, vceil 1e8, `dfloor_keep_velocity = true`.
  - Sponge: `he_sponge_rate 1e-3`, `he_sponge_r0 2.2265e13`, `mode radial`, `he_sponge_dmax 2e-14`. It acts on the
    whole floor plateau. `he_bc_outer = outflow`.
  - A DeltaAI twin (fl15, other floor/sponge) made the 1.03-1.10 R_ph layer worse. Decision 10-18: keep the
    baseline floor/sponge (Raven hg145_prod_1018/NOTE.md).

### 2.2 First steps of the 3-D run
- **Radiation start.** The same `ic_rad = fs_relax` with `ic_rad_relax_dtfac = 1`. The Raven gate log shows the
  ic_rad transient: L_top 5.8e4 Lsun at t = 50 s, settling to 1.04 L_MESA.
  - On one IC variant the first (cold) implicit solve did not converge at cycle 0, at the column/sph1d photosphere
    seam (job 31206670: dt_retry dt 100 -> 50 s, regrown by cycle 14). It did not recur.
- **Time stepping.** hesdirk2, `time2_start = be`, dt_retry on, `cfl_number 0.3`. dt is base-limited: 101 s on the
  wedge, 142.9 s on cs N=48. No dt or L ramp.
- **Scaffold.** `mlt_closure = adaptive` (frozen failed: L_top 0.39-0.58 L by 1 d). `mlt_flux_frozen = true`,
  explicit ramp `mlt_ramp_start 1.2856e6 s` (14.9 d), `mlt_ramp_time 6.4279e6 s` (74 d).
- **WB.** `wellbalance_dynamic = true`, `wb_option = polytropic`, `wb_cache_every 0` (rebuilt every stage),
  `wb_rmax 8.768e12` (1.1 R_ph).
- **Seed.** `he_seed 1e-3`, 16 modes, `he_seed_rad = true`, r 5.34e11-7.48e12. Wedge: k 1-2. cs: `he_seed_ylm =
  true`, l 1-8 (deck hgglob_prod_viper_ylm.athinput). HG262 uses the n2 band.
- **Background gates before production (the useful check for you).** Seed 0, full layout, 1-3 d.
  - Wedge viper bg chain (3 links, 1.27 d): interior quiet, 0.9-0.99 R_ph -0.37 km/s, 1.03-1.10 R_ph rho +32 %
    and slowing.
  - cs screen (Raven gpudev 31289336-38, 2 d): L_top/L_in 0.993-1.018 after 0.5 d. 1.1-3 R_ph v_r -0.05..-0.02
    km/s (sponge-terminal), max|dln rho| 2.8e-3. 0 retries, NC 0. VERDICT CLEAN.

### 2.3 Early trouble
- **Seam-driven photosphere restructuring.** With the wrong photospheric kappa and a frozen scaffold, the photosphere
  restructured (section 2.1). Fixed by the MESA-kap splice + adaptive closure. The 704 2-d screen was then clean:
  L_top/L_in 1.00 from 0.35 d, all bands <= 0.8 km/s.
- **The early HG expansion and the "ejecta" beyond R_ph (|v_r| up to 240 km/s, dM/M -8.5e-3) are a numerical
  start-up transient, not an outflow** (user 11-10). Exclude the first phase from statistics.
- HG did not need the floorcap (vceil 1e8, cool 2400 K floor, sponge on). I found no floorcap keys in the HG decks I
  read (uncertain whether a later HG run.cfg added them).

## 3. Old implicit-M1 handling of thin / floor gas (evidence from the decks and fork/prod-1114 source)
- There is **no explicit post-solve |F| <= cE clip** in the implicit path that I could find. `implicit_trans_limit =
  none` in both decks; lp/lpr exist but were not used. Realisability came from these pieces:
  1. The VET closure from the formal solution (`closure = vet_col` + vet_gd; f and g_top physically bounded).
  2. Face fluxes in `M1BeamFace` (rad_m1_implicit.cpp ~133-170). |f| is clipped to [0, 1], the M1 eigenvalue wave
     speeds are bounded by chat, and the HLL with a tau-dependent alpha uses the guard (1 - fbar^2). The blend
     weight `w_f` comes from `implicit_thin_*`: `implicit_thin_switch fse`, `_corr hc`, `_hcface tau`, `_js_face
     hcface`, `_recon jse`, `_bq lin`, h0/h1 0.3/0.6.
  3. `implicit_g0_limit = 1` (first-pass exchange clip to the energy a cell can exchange in one step),
     `implicit_g0_exchange`, `implicit_pos_gas`, `implicit_pos_floor` (positivity set, default since 10-02).
  4. `implicit_thin_freeze = 1e-2` (cells with c dt rho kappa_P < 1e-2 keep their opacities over the Picard passes).
  5. The floor-gas hydro guards: `implicit_cs2_clip = rad`, `implicit_wb_vkick_reject 1e9`, then the floorcap.
- The IC itself starts F/cE <= 1, because sph1d sets E from an RE transfer solution with the true F_r. fs_relax then
  relaxes E with F free but at the physical dt.

## Lessons relevant to your clean-scheme start
1. Do not start with a dt = 1e6 x dt_CFL BE pseudo-step in floor gas. Our working start relaxed the radiation at the
   run's own dt (`ic_rad_relax_dtfac 1`, at most 40 steps, hydro frozen, no radiative force on the gas, energy
   exchange on). Even so, L_top dipped to 0.75 L and recovered to 0.98 L.
2. Build E above the photosphere from a 1-D RE transfer solve with the column's own F (sph1d, taua 1 / taub 3,
   re_atm). F is then L/4 pi r^2 and E is consistent with it, so |F|/cE <= 1 at t = 0. The isothermal or
   constant-T guess for E is 2-2.6x off above R_ph (pgen comment, AG Car A).
3. The column's floor must equal dfloor (the R136 dfl17 change). The pgen's HSE/atmosphere marches stop at the
   column's minimum rho. A denser constant cap free-falls at (1-Gamma) g and drives the start.
4. Expect floor gas to heat and run to vceil in a no-line-force run. The fix that worked is the conservative floor
   cap (rho < 1e-14: |v| <= 300 km/s, T_gas <= T_rad, excess to E_rad). It restored a 3.5x larger dt and hesdirk2
   stepping, with a 1e-5 effect on the Fe-peak shell density. Without it, the vkick rejects forced 100 % BE.
5. Use wellbalance_dynamic with the background rebuilt every stage (wb_cache_every 0), not a reference frozen at
   t = 0. Put the WB reference only up to ~1.1 R_ph (wb_rmax), and keep `he_wb_grav_window` on.
6. Near-Eddington photospheres (Gamma 0.97-0.99) cannot be made static to 1 % (IC_HSE.md). Our runs accepted an O(0.1 g)
   seam residual and let the 3-D run inflate. The HG lesson: the opacity at the photosphere must match the source
   model, and the scaffold must hand over live (adaptive). Otherwise the photosphere restructures whatever the
   start.
7. Run a seed-free background gate (1-3 d, full layout) before production, and watch 1.0-1.1 R_ph drho and v_r,
   floor-band v_r, and L_top/L_in.

-- viper (for the USER)
