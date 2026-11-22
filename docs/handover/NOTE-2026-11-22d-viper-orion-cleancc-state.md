# NOTE 2026-11-22d viper -> orion (reply to TASK 2026-11-22a): clean-cc (cell layout) state for the momentum comparison

Viper has done no clean-cc work since the USER stop on 11-19 (NOTE 2026-11-20b c8d85bd3f has the full findings). Summary
and pointers:

## 1. Commit to run, and keys
- **Recommended: `fork/clean-cc-lbot-xtcons-viper` = 2eda97438** (pushed now, side branch). It is 65f075a38 (the lbot
  base: imposed-flux boundary-face stress fix) + our r9 work. With all our extra keys OFF it is the one that ran 1000 RSG
  restart cycles without collapse (run `lx_baselong`).
  - Exact command and keys: `/viper/ptmp2/jinma/ccfull_1117/runs/lx_baselong/cmd.txt`. Binary `ccint_1118/bin/he_star_m1_2eda97438`.
  - Key files: `/viper/ptmp2/jinma/ccfull_1117/keys/` (`rsg_prod_viper.txt`, `hg_prod.txt`, `hg_clean_fix.txt`,
    `r136_clean_*.txt`).
- **Also pushed:** `fork/clean-cc-lbot-pos-viper` = fdc30cf54. It is 2eda97438 + opt-in xt_cons, rst_k1, thin_low,
  thin_ydonor, low_face_d, bcg_fallback, the MPI deadlock fix (one MPI_MAX for the Picard exit) and diag_neg/ledger
  diagnostics.
  - **CORRECTION:** neither 2eda97438 nor fdc30cf54 has the Picard-exit MPI deadlock fix, so multi-rank runs can hang.
    Cherry-pick your own fix 3369185bc (CleanMaxRanks, in caltech-prod-bf4f917d4), or use the branch below.
  - **Also pushed now:** `fork/clean-cc-int-viper` = **06c9ee3cc**, a fast-forward of 6f3aa7d07 adding our deadlock fix.
    This int base DOES contain 65f075a38. It is what the HG restart below ran (as 6f3aa7d07).
- **Do not use:**
  - `fork/clean-cc` b53e35565 (push 18): it lacks 65f075a38.
  - The r9 base 245ea4418: L_bot comes out at -1e40 on RSG.
- There is no `docs/clean_*/{hg,rsg}_latest.txt` for cell. The cell keys are `rad_m1/scheme=clean clean_layout=cell` plus
  the lx_baselong cmd keys. Restart from an old prod rst needs **`clean_rst_k1=2`** (imports K1 without the face channels;
  =1 with faces gave ledger +0.2).

## 2. Thin transport, realizability, momentum (from the header of `src/rad_m1/rad_m1_clean_cell.hpp`)
- **Unknowns:** cell E and y = F/c. One implicit stage solves E and y together with the mixed-frame sources (Ma+2025 eqs.
  56-57), exact local exchange, gas T by bracketed Newton.
- **Face fluxes:** central average + HLL-type dissipation:
  - `fE = c(y_L+y_R)/2 - al (s/2)(E_R-E_L) + RC`
  - `fY = c(D E)_avg - al (s/2)(y_R-y_L)`
  - `s = c sqrt(max D_dd)`.
- **RC** is a Rhie-Chow-like momentum interpolation, `(1-al) c/chi_f [...]`, which turns fE into the compact diffusion
  flux where h c chi >> 1.
- **Dissipation weight al:** from the ray/tau sensor (`clean_al = tau | tau_sensor | tau_ray`). chi*dx enters only here,
  per the USER rule.
- **Reconstruction:** first order, or `clean_recon=plm` (van Leer weights lagged from the Newton iterate).
- **Thin transport = clean_consistency:**
  - Each face carries the bracket `C = F_ray - f(Ur)`, the rays' kinetic flux less the face flux of their own cell
    moments, so the ray state is a steady state of the discrete moments.
  - The momentum deviation `R = y - hE` streams with the rays' half-range factors `c(h+_L U_L + h-_R U_R)` and is
    absorbed at chi.
  - Rays: conservative FV short characteristics (`clean_rays = fv1|fv`) on the geodesic set.
- **Realizability:** no AP-HLL, no |F|<=cE clip, no FCT.
  - Positivity comes only from the stage-level BE fallback / dt retry when E<0 or a written-back eint<=0 appears.
  - `clean_thin_low=1` (faces with chi d<1 and hc chi>1 start low order) removed the first-BE E<0 cells. It did not
    stop the eint<=0 collapse, see section 3.
- **Momentum:**
  - Gas takes exactly the opposite momentum and total energy of the radiation source.
  - E changes only by face fluxes, so total energy telescopes to round-off.
  - Spherical polar: y holds physical components, and the rows use angular-momentum-conserving weights
    `r^-3 d_r(r^3 P_ra)` and `sin^-2 d_theta(sin^2 P_tp)`.
  - Cubed sphere: E and **L_z telescope** (the y_phi row is the torque about z, the conservation form of L_z = V R y_phi).
    y_r and y_theta are projected divergences, so linear radiation momentum does NOT telescope component-wise on cs,
    though a uniform isotropic stress is exact.
  - A cs energy/momentum balance must take the boundary fluxes from CleanCellRows (the solve's face areas): pcoord
    area x f0x1 is off by ~0.6 on cs.

## 3. RSG/HG restarts: the USER is right, they never ran cleanly
RSG runs: rsg.01658, 3 nodes x 2 MI300A; logs in `/viper/ptmp2/jinma/ccfull_1117/runs/<arm>`.

- **lx_baselong (2eda97438, extras off), 1000 cycles: stable but ~50x prod cost.**

  | | clean-cc | prod |
  |---|---|---|
  | dt mean | 94.5 | 931 |
  | s/cycle | 17.5 | 3.40 |
  | retries | 116 (112 NON-CONVERGED, 4 eint<=0) | — |
  | linear its/stage | 72 | — |
  | ledger | -9.7e-5 | 4.2e-5 |
  | L_top | 6.50e38 | 6.45e38 |
  | L_bot | 6.36e38 | 6.85e38 |

- **First BE step at full dt:** 510583 E<0 cells deep inside (rho 1e-7..1e-5). xt_cons leaves 425, in a thin shell at
  i~407-450. thin_low removes them.
- **xt_cons long runs (lx_long, lx_xtthin) collapse after ~100 cycles to dt 4-15 s.**
  - Almost every retry is a write-back eint<=0: the exact exchange has no root with eint>0, because transport drains
    more than the cell holds.
  - Ledger +0.06..0.10, L_top 20-28 % low.
- **Frozen opacity:** with `clean_opac_passes=2` the converged state depends on the first iterate. Every-pass
  opacities (1000) make xt_cons = base but give 21-cycle dt mean 20; untested long.
- **HG (`ccfull_1117/runs/s3_hgC0_i6`, binary 6f3aa7d07):** a 10-cycle restart took 1221 s. 9 of 10 steps were
  backward-Euler fallbacks (3 stage fallbacks), and dt fell to 0.90 s. It was never run longer.
- **CPU:** Raven CPU clean RSG/HG smokes timed out at 1.5 h inside cycle 1, because the cold start is too slow.
- **Restart pitfall:** old prod restarts carry 37 slope channels and clean keeps 136, so the first clean step is BE at
  full dt.
- **Keys:**
  - `time2_k1_redist` acts under scheme=clean, although the preflight says it has no effect.
  - `time2_face_fct` is forced off.
  - Local fix: rad_m1_clean.cpp:861 set dtr_fail without dtr_why ("reason ?"); dtr_why=4 is in fdc30cf54.

-- viper (for the USER)
