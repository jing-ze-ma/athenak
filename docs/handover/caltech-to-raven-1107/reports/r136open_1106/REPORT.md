## R136a1 sponge / open-outflow test: final report (11-06/07)

**Verdict: keep the production restart as it is.**
- No sponge setting raises dt: dt is set at 26-31 Rsun, below the sponge.
- No sponge (B) is slightly worse (-3 % dt) and the atmosphere falls in.
- he_sponge_mode=all (M) matches the control on dt and removes the vceil / table-top-T v_t cells.
- The infall without the sponge comes from a pgen bug, not physics. It needs code, so I stopped there (built a diagnostic binary, no merge).
- **REPORT.md was not written:** the tool refuses report files from subagents. The content is below; please save it as `r136open_1106/REPORT.md`.
- **Plot:** `/resnick/groups/carnegie_poc/jingze/r136open_1106/r136open.png`.

### Setup
- All arms restart from production rst 00005 (t=5.000 d), with production keys plus `rad_m1/time2_adm_poslim=true` and `bin_double=true`.
- Each ran 2 ranks on 1 GPU, about 3.1 s/cycle. That is not comparable with production's 1.15 s/cycle on 2 GPUs.
- Control is production's own hst/log from the same rst.
- **Noise member C:** perturbed rst (E×(1+1e-10u)), production binary and keys, run to 5.2 d. C vs production gives the noise: L_top rms 4e-4, Etot 4e-10, dt identical to within about 1 %. The hot v_t cells in C come and go (v_t 940-8200 km/s, table-top-T cells 1-12).
- Runs live in `/resnick/groups/carnegie_poc/jingze/r136open_1106/runs/`.

### Keys (pgen he_star_m1)
- **No sponge:** `problem/he_sponge_rate=0`. It returns before both the radial sponge and the dmax branch.
- **he_bc_outer=outflow (already in production):**
  - ghost rho ~ r^-2, outward v_r is copied, inward v_r is mirrored (not zeroed);
  - `he_wall_zero_flux` removes any inward mass, v_t and energy flux at the top face;
  - so the top lets mass out and blocks inflow.
- **A fully open top (inflow allowed):** needs code.
- **Arm M:** `problem/he_sponge_mode=all` damps v_t as well.

### Results, 5.00-5.15/5.20 d (no FATAL or non-converged steps in any arm)

| arm | mean dt vs C | dt dips | atmosphere | v_t / table-top-T cells |
|---|---|---|---|---|
| B (no sponge) | 5.47 vs 5.65 s (-3 %) | to 1.33 s | whole atmosphere falls at about 1300 cm/s2; top v_r -158 km/s, rms v_r 50-90 Rsun 158 (C: 22) | v_t reaches vceil, up to 20 table-top cells |
| M (mode all) | 5.66 s (= C) | 2.25 s once | held as in C | max v_t under 1200 km/s, 0-2 table-top cells |
| BL (no sponge + linear source, diagnostic binary) | -2.7 % | to 0.97 s | falls at about the physical rate (top -39 km/s at 5.2 d) | hot cells as in C |
| ML (mode all + linear source) | 5.61 s (= C) | | rms v_r 50-90 Rsun decays to 9 km/s | v_t under 900 km/s, 0-1 table-top cells |

- **Common to all arms:** R_ph, deep T and mass match C, and the E_tot budget closes.
- **L_top:** about -0.1 % at 5.2 d in BL and ML, which is 2.5x the noise.
- **Sponge power:** P_sponge about 4-6e33 erg/s, 2e-6 of L.

### Force imbalance above 40 Rsun: cause found
- **Size of the excess:** measured dv/dt with the sponge off is about -1300 to -1700 cm/s2 at every radius from 40 to 100 Rsun. Expected -g(1-Gamma) is -1342 at 40 Rsun and -260 at 99 Rsun. The excess grows with dr: zero where H/dr is about 1, about 1000 cm/s2 at H/dr 0.08.
- **Ruled out (1080-s tests from rst 5):**
  - WB reconstruction everywhere (`hydro/wb_rmax=1e15`, T1): same as T0 within 3 %.
  - Isothermal walk (T2b): same within 4 %.
  - Radiation kick off (`dbg_gas_force=false`, T3b): only 95-310 cm/s2 less inward, so the residual radiative force is small (a_ref is about 7 % above the code's κF/c).
- **Cause:** `HeStarGravity` (pgen line ~3251) applies the WB walked-pressure gravity source in every cell. It ignores `hydro/wb_rmax`, while the hydro reconstructs plainly above wb_rmax (27.5 Rsun). Where H is much smaller than dr, the walked face pressures over-pull the gas.
- **The sponge was hiding it:** the sponge was holding this force at about 1000 cm/s2. That is why the gas hangs and why it falls in blocks below the dmax switch.
- **Fix test:**
  - Patch in `/resnick/groups/carnegie_poc/jingze/r136open_1106/build/exp/src/pgen/he_star_m1.cpp`: opt-in `problem/he_wb_src_linear_rmin`, which uses the linear source -ρ(g - a_ref) above that radius.
  - Binary: `/resnick/groups/carnegie_poc/jingze/builds_rstfixx/bin/athena_gpu_he_star_m1_dc2dc20c_srclin_fmadoff`.
  - Without the key it is bitwise identical to production (T6, all hydro and M1 fields).
  - With the key (T5): -267 vs -260 expected at 99 Rsun, -379 vs -307 at 90.
  - Still open: 50-70 Rsun remains 1.4-1.9x too strong. Likely gas falling in from above or the opacity estimate; not tested.

### Production E_tot drop at 6.641 d (read-only check)
- **It is a real loss, not an hst artifact.** At hst row 573754→574000 s, E_rad fell by 6.77e42 erg and e_gas by 6.4e41. Mass, L_top and Mdot_top are unchanged, and the rows afterwards are smooth at the new level.
- **Size:** 6e-5 of E_rad, about 3300 s of L_top.
- **Trigger, cycle 131896:**
  - a `dt COLLAPSE` to 0.134 s on a floor cell: r = 35.4 Rsun, rho = dfloor 1e-17, p = 7.7e-3, printed T 5.2e14 K, cs 357 km/s, v = 0;
  - then 4 time2_posfix rounds that widened to 33,207 cells without the stage-1 increment;
  - then NOT ADMISSIBLE and the backward-Euler redo.
- **Other events:**
  - A smaller one of the same kind at 6.997 d: -4.8e41 in E_rad (posfix 37k cells, then 8.07e6 cells at cycle 137913).
  - The 8e6-cell posfix events at cycles 124584 and 124706 left no jump.
  - 7 dt collapses so far, all on hot low-density cells at 35-40 Rsun.
- **Not established:** where the energy goes (the BE redo vs the E-floor or rollback accounting). The end-of-run counters_rollback line should show it. Production was not touched.

### Decisions for the user
1. Restart choice: keep production as is, or switch to `he_sponge_mode=all` at a later restart (clean, same dt, hot v_t cells gone).
2. Whether to fix the WB gravity source above wb_rmax properly (honour wb_rmax in HeStarGravity). It needs GPU gates. The same code runs in HG and RSG wherever their tops have H much smaller than dr.
3. Whether to look into the posfix / BE energy loss on hot floor cells in production.

### State and cleanup
- **Allocations:** all arms are finished and all locks are released (02-03 GPU 0; 01-05 GPUs 0,1). Nothing of mine is running.
- **Kept:** final rst of B and M (5.15 d, 13 GB each, in `runs/B/rst` and `runs/M/rst`). The other test rsts are deleted. The run directory is 69 GB in total, on a filesystem that is 97 % full.
- **Left for you:** `runs/B_aborted1` is an 80 KB stub (my mistaken first B run when the target was shortened). Delete it if you like.
- **Scripts in `r136open_1106/`:** `an.py`, `series*.py`, `wbforce.py`, `force.py`, `tcmp.py`, `plot.py`, plus the queue scripts.
