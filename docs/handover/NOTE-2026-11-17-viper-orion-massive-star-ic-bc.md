# NOTE 2026-11-17 viper -> Orion (USER): how old production sets up a new massive-star run (answer to 16b124e28)

This note answers QUESTION-2026-11-17-orion-viper-massive-star-ic-bc.md. I did not work from memory. Every key below was read
from one of these sources:

- **Restart input blocks** (the `<par>` text at the head of each file): R136
  `/viper/ptmp2/jinma/r136fresh_1108/run/rst/r136.00096.rst`, HG `/viper/ptmp2/jinma/prod_1107/hg/rst/hegiant.00171.rst`,
  RSG `/viper/ptmp2/jinma/prod_1107/rsg/rst/rsg.01545.rst`. These files were only read.
- **Fresh-start decks:**
  - R136 arm E2: `/viper/ptmp2/jinma/r136fresh_1108/decks/r136direct3d_e2_viper.athinput`;
  - M35MMS: `/viper/ptmp2/jinma/m35mms_1115/decks/m35mms_ours_mesaR.athinput`.
- **run.cfg KEYS:** `prod_1107/{hg,rsg}/run.cfg` and `r136fresh_1108/run/run.cfg`.
- **IC file headers** (provenance and the ic_star spec JSON).
- **Run logs.**
- **Source:** `fork/prod-1118` (b2d852bb6, which contains the 1d1972805 production binary), `src/pgen/he_star_m1.cpp` (HS)
  and `src/rad_m1/rad_m1.cpp`.
- **Worker notes:**
  - `r136fresh_1108/{DECKS,IC_HSE,RELAX}.md` and, on Raven, `/raven/ptmp/jinma/r136fresh_1108/{IC_HSE2,DTCOLLAPSE}.md`;
  - `m35mms_1115/SETUP.md`.

None of these are restarts used for a test, so no rotation applies. All four runs are fresh t = 0 starts from an IC column.

## TL;DR

1. **The relax includes MLT.** The he_ic_equil relax runs with the frozen MLT scaffold at w = 1. It is deposited as
   `<rad_m1>/esrc`, so radiation carries L_in - L_MLT, and the force balanced by the march is the force of that relaxed
   radiation.
   - In old production the residual inside the scaffold is 1e-5..1e-4 g. For R136, after 1 cycle it is max 4.5e-5 g below 24.5 Rsun.
   - Your 0.1-0.4 g in cells 50-150 is therefore what we would expect if your relax leaves out F_MLT.
2. **The Gamma > 1 layer is not balanced. Nobody balances it.**
   - In R136 it sits at 25.04-25.31 Rsun (1.742e12-1.761e12 cm), which is exactly your ~2 g at r ~ 1.74-1.76e12.
   - Old production has the same residual there: 2.1 g at active cell 155 (25.04 Rsun) in round 1. The march stops
     there, and the cells above keep the pgen state.
   - More rounds diverge: maxit 4 gave 2.1 -> 0.8 -> 1.0 g.
   - IC_HSE2 shows that no static column exists there with the code's opacity, L = L_in and the IC's scaffold. Production
     accepts the residual and lets the layer inflate.
3. **Production uses maxit 1 with last_march = true in all four runs.** The relax is 40 steps at 1 dt
   (he_ic_equil_dtfac = 1.0, relax_tol 1e-4). There is no separate radiation-only phase, no dt ramp and no damping. The only
   time-dependent start-up device is the auto MLT scaffold ramp, w(t) cosine 1 -> 0 over t_sat.

## 1. IC source and mapping

**Tool and mode.** Every run uses the same tool: `tools/ic_star make_star_ic.py`, mode `mesa`, scheme `mesa`, 7-column file
`r rho eint F_r(radiative part) E | T F_MLT/F`, read with `he_ic_cols = 5`. The doc string of `icstd/mesacol.py` describes the
mode. A copy is at `/viper/ptmp2/jinma/m35mms_1115/tool_ic/icstd/mesacol.py`.

**Structure.**
- T(P_tot) follows MESA's ln T(ln P) through a PCHIP fit, so the temperature gradient nabla is MESA's own, including MESA's
  MLT/MLT++ where MESA used it.
- Hydrostatic equilibrium is re-integrated OUTWARD from the column bottom with our EOS: P_gas,code(rho, T) = P - aT^4/3.
- Gravity is the code's: G m(r) = G m(r_bot) + G int 4 pi r^2 rho dr. The deck has `he_gm` = G m(x1min) and
  `he_gm_column = true`.

**Flux split.**
- F_rad/F = nabla_MESA / nabla_rad(our kappa, g = G m/r^2, F = L/4 pi r^2).
- F_MLT/F = 1 - F_rad/F, clipped to [0, 1].
- Where F_rad > F, the column takes F_r = F, the excess is reported, and those layers lose heat in the run.
- So the scaffold is the residual of MESA's structure. It is not an MLT formula of ours. Spec `use_mltpp: false`; the
  Cox MLT keys alpha 2 (HG, M35) and 3 (RSG) are used only outside the MESA range.

**Top of the column.**
- The column stops where T(P) reaches the Hopf tau_j = 3 temperature.
- A Hopf layer is then integrated outward to tau 2/3, with gas HSE and the radiation force of F. This gives our R_ph and
  rho_ph.
- Above R_ph comes an isothermal gravity-only atmosphere at T_ph, down to the column floor, then the floor up to the column
  top.
- The column floor rule is `rho_floor_rel 4e-6` x rho_ph, rounded to a decade:

| run | column floor (from the IC file) |
| --- | --- |
| HG | 1e-14 |
| RSG | 1e-17 |
| M35 | 1e-15 |
| R136 | 1e-17 (after `mkic_dfl.py`, see the R136 row below) |

**Seam pass.**
- HG and RSG: on (spec `seam.enable true`, maxit 8). The CPU code is run on the column and its dump is fed back.
- R136: the IC header reads "seam fixed point (seam_only.py), seam pass 2: ... scaffold of the dump (taper off)".
- M35: off (SETUP.md sect. 9).

**Per-run IC inputs:**

| run | MESA model | IC file (md5 where known) | opacity | R_ph, Teff (IC header) |
| --- | --- | --- | --- | --- |
| R136 (sp) | Keszthelyi+2025 TS MLT++ compact, 233.9 Msun, log L 6.833 (Caltech r136a1_1104) | `r136fresh_1108/pkg/ic/ic_r136a1_B_dfl17.txt` (e8de8c76). It is Caltech's ec3c8a81 with the 2517 cap rows lowered from rho 1e-13 to 1e-17 = dfloor (`jobs/mkic_dfl.py`), so the code's HSE atmosphere runs down to dfloor instead of ending in a 1e-13 cap that free-falls | Rosseland `rosseland_mesakap_a09_x0.58_z0.0064_topsbase_T7`, Planck `planck_tops_..._T7` | 24.91 Rsun (MESA), 59.1 kK |
| HG (cs) | S. Wu 14p5M_widebinary_profile14, He giant | `prod_1107/tables/hgglob_1021/files/hg145-prod-1018/ic_hegiant_14p5M_mk.txt` | `rosseland_mesakap_afe2dd8da2bd` (mesa_kap join log T [4.6, 4.8]), `planck_ext2_..._416a6119` | 114.2 Rsun, 7.50 kK |
| RSG (cs) | `rsg_g25_1019/mesa/m22_p438_sel.data` (anchor: mass) | `rsg_E_1020/files/rsg-E-1020/ic_rsg_g25_deep.txt` | `rosseland_mesakap_4c0e8b64d26a` (join [6.3, 6.5]), `planck_ext2_..._d38e5526` | 833.7 Rsun, 4.11 kK |
| M35MMS (sp) | Will Schultz `35M_save_profile_MS.data` | `m35mms_1115/ic_ours_mesaR/ic_m35mms.txt` (6e554487) | MESA kap over the full range (join [8.7, 8.8], `rosseland_mesakap_d0757106cce6`) + TOPS ext2 Planck (`59064cc8e7ba`). This is the USER 11-15 default for new runs | 13.69 Rsun, 35.8 kK |

**On the grid** (HS ~650-700 and he_ic_balance ~1883):
- The column is log-interpolated on `he_nfine` fine nodes. M35 uses 262144; the default is 16384.
- With `he_ic_eint_from_t = true`, eint = e(rho, T_col).
- `he_ic_rad = sph1d`, `he_ic_rad_hse = true`, `he_ic_rad_gas = re_atm`, `he_ic_rad_atm_force = module`,
  `he_ic_rad_taua = 1`, `he_ic_rad_taub = 3`: a radiative-equilibrium re-solve above tau_R = 1 using the module's own force;
  below tau 3 the column T is kept.
- `he_ic_rad_hse_gmax = 0.9` (default): the pgen HSE march clamps the effective gravity at 0.1 g wherever Gamma > 0.9.

**Hydro floors (`<hydro>` dfloor, tfloor_kelvin):**

| run | dfloor | tfloor_kelvin |
| --- | --- | --- |
| R136 | 1e-17 | 16500 |
| HG | 1e-18 | 2400 |
| RSG | 1e-20 | 2000 |
| M35 | 1e-19 | 24800 |

**Velocity and seed.** The initial velocity is zero. The only perturbation is `he_seed`, applied after the equilibrium
(the deferred seed):
- R136 E2: `he_seed = 1.0e-2`, `he_seed_band = old`, `he_seed_kmin 2`, `he_seed_kmax 12`, `he_seed_ylm = false`,
  `he_seed_rlo/rhi = 1.4345e12/1.8436e12`, `he_seed_rng 1104`.
- HG, RSG and M35: `he_seed = 1e-3`, n2 band (rlo/rhi in the decks).

## 2. he_ic_equil in practice

**Keys.** All four runs use the same keys:

    he_ic_balance = true   he_ic_equil = true   he_ic_equil_maxit = 1   he_ic_equil_last_march = true
    he_ic_equil_dtfac = 1.0   he_ic_equil_relax_n = 40   he_ic_equil_relax_tol = 1.0e-4   he_ic_equil_tol = 1.0e-9
    <rad_m1> ic_rad = fs_relax  ic_rad_fs_iter = 3  ic_rad_relax_n = 40  ic_rad_relax_tol = 1e-4  ic_rad_relax_dtfac = 1
    mlt_flux_frozen = true   mlt_closure = adaptive   mlt_closure_taper = 0.01

The code defaults are different (maxit 30, dtfac 1e6, relax_tol 1e-12); production names the values above.

**What one round does.** From the HS comment at the enrolment (~2968), checked against `HeStarIcEquil` (~5123):
1. `ic_rad` relaxes E, F and the gas T at fixed rho. The relax includes the frozen MLT scaffold at w = 1, i.e.
   **L_rad + L_MLT = L_in at every face**, and then makes one pass with the radiative force on (`icr_fmeas`).
2. a_ref = f/rho from that deposit. Phi_eff is rebuilt, and rho is re-marched at the relaxed T in the discrete WB face
   balance, OUTWARD from the wall cell (whose rho is kept). Inward marching is unstable.

**How the scaffold enters.** `mlt_flux_frozen` deposits F_MLT = F_r fmlt/(1 - fmlt) (column 7) as the conservative source
-(A F_MLT|+ - A F_MLT|-)/V into the RADIATION energy through `<rad_m1>/esrc` (HS ~660). So yes, the march balances gravity
against the force of L_in - L_MLT.

**Where the march starts.** The he_ic_balance anchor is `he_ic_balance_rmax`; cells above it keep the column state:

| run | he_ic_balance_rmax |
| --- | --- |
| R136 | 1.7178e12 = 24.69 Rsun |
| HG | 7.0208e12 = 100.9 Rsun |
| RSG | 5.1193e13 = 735.8 Rsun |
| M35 | 8.389e11 = 0.8808 R_ph ("inside the MLT++ adiabat, below the tau-3 join") |

**Measured round-1 numbers (from the logs):**
- **M35** (`prod_viper_ours_mesaR/run.log`, first start):
  - module residual max |a_net|/g = 0.112 at r = 9.579e11 = 13.77 Rsun, just above R_ph 13.69;
  - WB face residual 1.27e-2; max |L(face)/L_in - 1| = 4.2e-3;
  - re-marched cells 0..265 of 336; max |drho/rho| 0.28 at 9.67e11 (above R_ph);
  - per-cell table in `he_ic_equil.txt` in that directory.
- **R136** (IC_HSE.md, Raven 32048340 arm A1, dfl17 IC). Startup residual R = a_r/g after 1 cycle:

| r [Rsun] | R |
| --- | --- |
| < 24.5 | max 4.5e-5 |
| 24.72-24.94 | +5e-4 .. -8e-3 |
| 25.03-25.25 | +0.13 .. +0.63 (one shell +1.9) |
| R_ph 25.37 | +0.25 .. +0.005 |
| floor 26.15-100 | -0.28 .. -0.33 = -(1 - Gamma_es) |

  - At active cell 155 (25.04 Rsun) the re-march finds no face solution: the round-1 module residual is 2.1 g, and
    |da/a| = 2.45 after 200 passes. The march stops there.

**The Gamma > 1 / iron-bump layer.** Production leaves it unbalanced on purpose. IC_HSE2.md (Raven) found:
- The scaffold ends abruptly at 25.04 Rsun: column 7 drops 0.667 -> 0 in one row. This comes from Caltech's seam pass
  ("taper off").
- Above that, all of L is radiative and kappa_R = 0.46-0.66 > kappa_Edd = 0.449.
- Measured Gamma_code = (f/rho)/g:

| r [Rsun] | Gamma_code |
| --- | --- |
| 25.05 | 1.33 |
| 25.12 | 1.17 |
| 25.22 | 1.06 |
| 25.29 | 1.006 |
| 25.32 | 0.999 |

- Below the top, inside the scaffold, Gamma_code is 0.53-0.60, with F_MLT/F 0.67-0.81 at 24.87-25.02 Rsun.
- Outward ODE tests with the code's opacity: from the code state the density inverts and Gamma runs away to > 5 within
  0.01-0.02 Rsun.
- An Eddington-capped frozen scaffold (eps 0.02-0.10) still collapses T by 25.11-25.14 Rsun. A frozen flux deposited into
  E is not a convective temperature gradient, so it cannot hold the layer up.
- Verdict there: a static IC does not exist. The recommendation, which production took, was to accept the inflation as
  physics and launch from the dfl17 IC.

**The other runs.**
- The IC is not modified to remove the layer anywhere.
- M35's IC has its Fe CZ inside the scaffold: 12.74-13.26 Rsun, max F_MLT/F 0.158.
- I did not check whether HG and RSG have a Gamma_rad > 1 layer above their scaffold top.

**At the top.**
- Cells above `he_ic_balance_rmax`, and above the cell where the march stops, keep the pgen state: the module-force HSE
  atmosphere, gmax-clamped at Gamma > 0.9.
- Floor gas cannot be static and falls at (1 - Gamma) g. In R136 that is -0.28..-0.33 g above 26.15 Rsun.
- Gamma clamp effect (R136): +(Gamma - 0.9) g, i.e. R = +0.064 at R_ph.
- Arms that were worse (IC_HSE.md): gmax 0.99 puffs the atmosphere to 26.6 Rsun; gmax 0.999 FATALs; maxit 4 and maxit 2
  diverge.

**After t = 0: the adaptive closure.** `mlt_closure = adaptive` (HS ~228):
- F_sub relaxes toward s(r) max(0, L - L_rad - L_conv,res)/A over `mlt_relax_time` (1500 s; RSG 1.5e4).
- s(r) = clamp(fmlt_IC/0.01, 0, 1), i.e. only inside the IC's FeCZ mask.
- F_sub(0) = the IC's F_MLT.
- w(t) caps it as an overall factor, and the closure is skipped once w = 0 (HS ~3644).

## 3. Boundaries

**Inner x1, radiation.**
- `implicit_bc_x1min = flux` with `implicit_flux_x1min` = the column F_r at x1min.
- he_star_m1 FATALs if the key differs from the column's F_r by more than 1e-3.
- Where the column carries F_MLT at the wall, the key is not L/(4 pi r^2). The ratio key/(L/4 pi r_in^2) per run:

| run | x1min | key/(L/4 pi r_in^2) |
| --- | --- | --- |
| R136 | 19.45 Rsun, log T 6.50 | 0.99942 |
| HG | 8.0 Rsun | 0.98763 (wall F_MLT 1.2 %) |
| RSG | 25.31 Rsun | 0.99740 |
| M35 | 9.697 Rsun | 0.99942 |

**Inner x1, hydro.** All runs: `he_bc_inner = wall`, `he_wall_closed = true`, `he_wall_zero_flux = true`.
- The ghosts are the (balanced) column scaled to the edge cell.
- `wall_closed_ix1`: the wall face sees the exact mirror state, so the mass flux is zero to round-off.
- `he_wall_zero_flux`: the inner face also loses its mass, transverse-momentum and energy fluxes; the pressure flux stays.
- Gravity: `he_gravity = wb` + `he_gm_column = true` (static enclosed-mass gravity of the IC column).
- `<hydro>`: `wellbalance_dynamic = true`, `wb_x1 = true`, `wb_option = polytropic`; `etotgrav` is required. Values of `wb_rmax`:

| run | wb_rmax |
| --- | --- |
| R136 | 1.916e12 = 27.5 Rsun |
| HG | 8.768e12 = 1.10 R_ph |
| RSG | 6.393e13 = 1.10 R_ph |
| M35 | 1.0477e12 = 1.1 R_ph |

- R136 also has `problem/he_wb_grav_window = true`. It is the code default since prod-1108 unless recipe = none.

**Outer x1, radiation.** All runs: `implicit_bc_x1max = marshak`, `implicit_flux_x1max = 0.0` (no incoming flux).

**Outer x1, hydro, and sponge.** The sponge (HS ~3855) relaxes v at `he_sponge_rate` x w, with w = ((r - r0)/(r_top - r0))^2
above r0, and w = 1 anywhere rho < `he_sponge_dmax`. With mode `radial` only v_r is damped. The removed KE leaves E_tot.

| run | he_bc_outer | he_sponge_rate [1/s] | r0 | dmax | mode |
| --- | --- | --- | --- | --- | --- |
| R136 E2 | `noinflow` | `0` (sponge OFF; coordinator 11-08: the density-threshold sponge made radial streaks) | (6.2613e12, inert) | (2e-13, inert) | radial |
| HG | `outflow` | `1.0e-03` | 2.2265e13 = 0.9 x1max | 2.0e-14 = 2 x column floor | radial |
| RSG | `outflow` | `1.0e-04` | 1.585e14 = 0.9 x1max | 2.0e-17 | radial |
| M35 | `outflow` | `1.0e-01` | 9.577e11 = 0.9 x1max | 2.0e-15 = 2 x column floor | radial |

- `outflow`: the ghosts continue the edge cell as a constant-velocity wind with rho ~ r^-2, and no inflow.
- `noinflow` is the code default.

**Lateral boundaries.**
- R136: sp wedge, theta pi/3..2pi/3, phi 0..pi/3, periodic x2 and x3 with `sp_x2_periodic_image = true`, 128 x 128.
- M35: Schultz+22 wedge, theta [4pi/9, 5pi/9], phi [0, pi/9], 336 log r x 256 x 256, periodic, image on.
- HG and RSG: cubed sphere, full sphere, N = 128 per face (`ix2_bc = panel`), 24 blocks 704/548 x 64 x 64. cs-only keys in run.cfg:
  - `mesh/cs_seam_blockend=extend`, `mesh/cs_vertex_3cell=true`, `mesh/cs_seam_rho_guard=4`;
  - `rad_m1/cs_cross_mono=true`, `rad_m1/cs_vertex_mono=true`, `rad_m1/implicit_cs_e0_metric=true`,
    `rad_m1/implicit_cs2_wbsolve=true`;
  - HG: `rad_m1/implicit_am_cs_kecap=0.5`, `implicit_clip_kecap=0.5`, `time2_inc_kecap=0.5`.

**WB and force reference.** `<rad_m1>` `force_reference = wb_arad`, `force_reference_work = split`. he_star_m1 FATALs
without wb_arad when he_gravity = wb.
- The WB pair carries rho a_ref of the equilibrium, so the coupling hands the gas only the residual radiative force.
- With split, the gas receives the reference work and E pays it ("exact" pay by default).
- a_ref is the he_ic_equil a_ref and stays fixed for the run. It travels in the restart as trailer `HSDICE01`; a restart
  FATALs if the trailer is missing or the grid differs.

## 4. Start-up procedure and transient

**In Driver::Initialize.** pgen column -> he_ic_balance -> ic_rad (fs_relax, 3 iterations) -> he_ic_equil round 1 with
the last march -> deferred seed -> cycle 1.
- That is all of it. There is no radiation-only settling phase, no dt ramp and no damping time. dt simply starts at the
  CFL value: `time/cfl_number = 0.3`; `cfl_rad 0.4` in M35.
- R136 E2 started at dt 15.4 s (DTCOLLAPSE.md arm E smokes, t = 153.85 s at cycle 10).

**The only time-dependent device: the MLT scaffold ramp.**
- `mlt_ramp_start = auto`, `mlt_ramp_time = auto` (n2ramp-1022).
- On a fresh start: start 0, time t_sat = 10.5/(0.61 sigma_max) from the IC column's N^2, stored in the restart (`HSRAMP01`).
- Over that time w(t) goes 1 -> 0 as a cosine. It multiplies both the frozen deposit and the adaptive F_sub, so after
  t_sat the scaffold is OFF and convection must carry the flux.
- Values printed in the run logs:

| run | t_sat [s] | t_sat [d] | log |
| --- | --- | --- | --- |
| R136 | 10312.3 | 0.119 | `r136fresh_1108/run/run.log` |
| HG | 1.15275e6 | 13.34 | `prod_1107/hg/run.1117_050458.log` |
| RSG | 9.09149e6 | 105.2 | `prod_1107/rsg/run.log` |
| M35 | 30544.8 | 0.354 | `m35mms_1115/prod_viper_ours_mesaR/run.1116_022925.log` |

- RELAX.md shows why the ramp matters. An 8 x 8 relax run kept w = 1 with the adaptive closure. The scaffold
  double-counted the resolved FeCZ convection: L_MLT 7e39 -> 3.3e40, L_rad 4 L_in, Mdot 0.4 Msun/yr.

**What "settled" means.**
- There is no settled-or-not gate before production. The initial transient is simply excluded from statistics.
- USER 11-10 on HG: the first expansion is a numerical transient, not an outflow. Science uses the relaxed envelope,
  checked for drift.
- HG relaxed values (memory hg-no-outflow-1110, from the HG page): R_ph 0.90 R_IC, T_eff 7.84 kK, L_top/L_in 0.985.
- R136 is not settled after 37 d. `r136expand_1116/ana/a3d_report.txt`: R_ph (shell-mean tau_R = 2/3) went 34.2 -> 76.4
  Rsun (window means 0-6 d and 30-37 d), slope 11 km/s overall.
- Start-up history of R136 (the alternatives that were tried):
  - The 1000-cycle start (IC_HSE.md, Raven 32048341) moved the photosphere outward. L_top fell 2.12e39 -> 0.74e39 while
    L_bot stayed 2.215e39.
  - Relax-then-remap (RELAX.md, `he_ic_remap_file`) was built, but arm D inherited floor pockets and was stopped.
  - Production is arm E2: a direct 3-D start from the dfl17 IC.

**R136 lateral convection finding** (`/viper/ptmp2/jinma/r136expand_1116`). The finding comes from the coordinator's
summary of that directory's arms. It has no written note there, and I did not re-derive the number:
- The 1-D/narrow-wedge tests show lateral convection at 20-25 Rsun with an e-fold time of about 0.09 d.
- The wedge erupts unless lateral motions are damped.
- Arms: `runs2/{w8,lb,lc,n05,n025,n0125,t2d_*,t2e_*,dir24x8,dir48x8,d48x8b}`, decks `decks2/`, a latdamp build in
  `build_lat/athena_he_1d197280_latdamp`.
- For your 1-D/half-phi clean runs, this means a balanced IC will not stay quiet. The FeCZ goes unstable on ~0.1 d.

## 5. Differences between the runs

| item | R136 (sp) | HG (cs) | RSG (cs) | M35MMS (sp, newest fresh start) |
| --- | --- | --- | --- | --- |
| geometry | wedge pi/3 x pi/3, 128^2, 831 r (deck). Since 12.15 d: g150w stretched grid, 738 r to 150 Rsun, IC `..._dfl17_ext165` | cs full sphere N128, 704 r, 8-355.6 Rsun | cs full sphere N128, 548 r, 25.3-2531 Rsun | 20 x 20 deg wedge 256^2, 336 log r, 9.7-15.3 Rsun |
| outer hydro | noinflow, sponge off | outflow + sponge 1e-3 | outflow + sponge 1e-4 | outflow + sponge 0.1 |
| scaffold ramp t_sat | 0.12 d | 13.3 d | 105 d | 0.35 d |
| mlt_relax_time | 1500 | 1500 | 1.5e4 | 1500 |
| IC seam pass | Caltech seam pass 2, scaffold "taper off" (abrupt top) | on | on | off |
| Gamma_rad > 1 above the scaffold | yes, 25.04-25.31 Rsun (IC_HSE2) | not checked | not checked | Fe CZ inside the scaffold; max residual 0.11 g just above R_ph |
| seed | 1e-2, k 2-12, band old | 1e-3 n2 band, Ylm, rng 22102 | 1e-3 n2 band, rng 1234 | 1e-3 n2 band (the sp run ignores he_seed_ylm) |
| extra keys | `rad_m1/floorcap_rho=1e-14`, `floorcap_v=3e7`, `implicit_wb_vkick_reject=1e9`, `hydro/eos_logt_max=7.0` | cs keys (sect. 3), kecaps 0.5, `eos_logt_max=8.0` | cs keys, `hydro/merge_lat_blockend=exchange`, `eos_logt_max=8.0` | HG262 key set (`KEYS_P1118_HG`) + S22 wedge |
| binary now | `athena_he_gpu72_x_1d197280_p1118` | same | rerun on Raven (pHG262 holder); I did not check its binary here | `athena_he_gpu72_x_0940c7d1_t9fix2` |

The radiation numerics keys (vet_gd_*, time2_*, implicit_*) also differ between the HG/RSG run.cfg and the R136 run.cfg.
They do not touch the IC or the boundaries; the full lists are in those run.cfg files.

## Unknown or disagreeing

- **RSG ramp key comment.** The RSG and R136 decks carry the comment "adaptive, no ramp" (rsgg25-1019) on
  `mlt_closure`, but the restart trailers hold auto ramp numbers (105 d and 0.12 d), and the logs apply them
  (w = 0 now). The comment is stale; the code and logs are what ran.
- **Orion's Gamma 1.69.** Your Gamma = 1.69 at 1.74-1.76e12 assumes all of L_in is radiative. Our measured Gamma_code
  peak there is 1.33 (tab 1.48) at 25.05 Rsun, with the dfl17 IC and the scaffold zero above 25.04. The difference may
  come from your relaxed T/rho; I have not checked it.
- **Your cells 50-150.** I do not know which grid you use (831 or 738). On the 831 grid, cell 155 is 25.04 Rsun, so
  cells 50-150 lie inside the FeCZ scaffold (F_MLT/F up to 0.81). That is consistent with the missing L_MLT in your
  relax.
