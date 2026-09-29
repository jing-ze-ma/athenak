# TASK for DeltaAI: WASP-121b 1x rot 300 -> resistive MHD, max_eta scan (bbot 3 G and 10 G)

**User 09-29:** continue the WASP-121b 1x hydro production from rot 300 as **resistive MHD at the same
resolution** (C32, nx1 76), with a bottom field of 3 G and of 10 G, and find an acceptable Ohmic cap
`<mhd>/max_eta`. The scan runs on DeltaAI (GH200). viper built the recipe and ran one smoke only.

## Package (made on viper; the user copies it by hand)
- Tarball: viper `/viper/ptmp2/jinma/w121_mhd_0929/w121_mhd_etamax_pkg.tgz`, 128 MB, md5
  12f574c79f1a056ffbae2d918f15fe0d. Put it on DeltaAI at `<PKGROOT>` (the user fills this in; project
  space, not home), then run `tar xzf w121_mhd_etamax_pkg.tgz`. This gives `<PKGROOT>/deltaai_pkg/`
  (171 MB, file md5s in `MD5SUMS`).
- Contents:
  - `dhj.00600.rst`: the viper w121prod_0929/w1x restart at **rot 300** (t = 3.304605e7 s, P_rot 110153.5 s,
    141 MB). It is the source state. It is not needed to run, because remap.dat is made from it.
  - `remap.dat`: the IDENTITY remap of that restart (`docs/handover/scripts/dhj_remap.py`, no `--grid`; the
    budget printout in `remap_identity_printout.txt` shows 0 change in every column, 0 cells at dfloor).
  - `remap_mhd.athinput.template`: the embedded input of the restart with these changes, and nothing else:
    `<hydro>` -> `<mhd>`; rsolver lhlld; `ohmic_resistivity = eos`, `max_eta = @MAX_ETA@`,
    `use_rkg_sts = @STS@`, `cs_lowbeta_fallback = 0.5` (as cs_mhd_prod4); `problem/bbot = @BBOT@`;
    `problem/ck_impl_conserve = 1`; output3 = `mhd_w_bcc` 4 per rotation; `remap_file`, `ic_profile`
    and ck tables point into the package (`exo_fms_ck/` and `ic_w121_1x.txt` are included).
  - `make_arm.sh ARMDIR BBOT_GAUSS MAX_ETA STS [PKG] [DFLOOR]`: fills the template and adds a `mhd_divb`
    output once per rotation.
  - `scan.sh RUNROOT [NROT] [submit]`: sets up the 16 arms below; with `submit` it runs one sbatch per arm.
  - `run_arm.sub`: the DeltaAI job, one arm on 1 node x 4 GH200. **Edit `@ACCOUNT@` and `@BIN@`.** The first
    link is a fresh start (`-i mhd.athinput`); later links restart from the newest `rst/` file. tlim is
    rot 300 + NROT, so all arms stop at equal physical time.
  - `ana_eta.py ARMDIR ...` (+ `bin_convert.py`, and `coarse.bin` = x_e/T table of the run composition):
    the analysis described below.
  - `viper_smoke/`: the viper smoke inputs, logs and analysis.

## Build
rt-integration >= fe7301e5, dhj target (PROBLEM = deep_hot_jupiter_rt). Use your incremental worktree
script (`build_inc_deltaai.sh`) with CUDA, Kokkos_ARCH_HOPPER90 + Kokkos_ARCH_ARMV9_GRACE and cray-mpich
GPU-aware, as for w121prod. MHD needs no build flag: `<mhd>` in the input selects it. Record the md5.
Launch as in NOTE-2026-09-28-deltaai-w121prod.md: all GPUs visible,
`KOKKOS_MAP_DEVICE_ID_BY=mpi_rank`, `SLURM_CPU_BIND=cores`, and no `--gpus-per-task/--gpu-bind`.

## Ranks
There are 24 MeshBlocks (6 panels x 2 x 2 of 16^2 x 76). Valid rank counts are 1, 2, 3, 4, 6, 8, 12 and 24.
**1 node = 4 GH200 = 4 ranks (6 MeshBlocks each) per arm** is the default in run_arm.sub. Two arms can
also share a node at 2 ranks each.

## viper smoke (job 12027359, apudev, binary rt-integration ae767d20 ROCm 7.2, 2 ranks, 50 cycles)
- Setup: bbot 3 G, max_eta 1e13, dfloor 1e-16 (production). There were three variants: (a) no STS,
  (b) max_eta 1 (eta ~ 0), (c) STS.
- All three: rc 0, 0 FATAL, 0 NaN, 0 NOT-CONVERGED. The remap file was read.
- Initial B: max |sum area*B|/(|B| max area) = 9.7e-15, and max |B| = 0.8456 code = bbot. After 50
  cycles, max |div B| dr/|B| = 9.7e-15.
- **dt = 0.456 s at cycle 0, rising to 0.560 s at cycle 50. The hydro dt is 13.8 s.** dt is identical in
  (a), (b) and (c), so it is **not Ohmic**: it is the MHD CFL. The limiting cell was identified from the
  cycle-0 bin, and the estimate reproduces 0.456 s. It is the density-FLOORED top: rho = dfloor = 1e-16,
  p = 3.6e-11 bar, where the dipole |B| = 0.45 code gives v_A = 450 km/s (dr = 700 km).
  - Estimated dt at 10 G: 0.14 s.
  - Estimated dt with a higher `<mhd>/dfloor`, same state (3 G / 10 G):
    - 1e-15: 1.3 / 0.39 s;
    - 1e-14: 3.2 / 1.1 s;
    - 1e-13: 7.4 / 3.0 s (7 % of the cells are raised, all at p < 3e-8 bar; added mass 4e-10 of the total);
    - 1e-12: 12.7 / 7.1 s (raised cells reach 2.7e-7 bar).
- **USER DECISION:** the main scan uses **dfloor 1e-13**. By the standing rule, the top (p < 1e-6 bar)
  does not count for accuracy, but it must not set dt. Two check arms keep the production dfloor 1e-16
  (b3_e13_f16 and b10_e13_f16).
  - Without a raised floor, every arm is v_A-limited at 0.14-0.46 s, and the Ohmic cap then costs nothing
    up to about 5e13.
  - A targeted alternative needs code: a magnetic density floor rho >= B^2/v_A,max^2 in the top. It does not
    exist; ask the user before writing it.
- Ohmic dt on the rot-300 state, as NewTimeStepGeneralResist x cfl 0.3:
  - 17.5 s at 1e12, 1.78 s at 1e13 and 0.185 s at 1e14;
  - the binding cell is at about 1.8 bar on the night side.
  - So 1e12 is free, 1e13 costs about 4x against the 7.4 s CFL (3 G, dfloor 1e-13), and 1e14 costs
    about 40x without STS.
- ms/cycle on 2 MI300A: 30.7 (MHD, including startup), against about 17.5 for the hydro.

## Scan (scan.sh): run 1-2 rotations per arm first (NROT = 2)
The main scan uses dfloor 1e-13, bbot 3 and 10 G. Arm names are b<G>_e<log eta>[_sts]_f<-log dfloor>.
- b{3,10}_e0_f13: eta ~ 0 (max_eta 1), ideal MHD with numerical resistivity only. This is the
  "cap -> 0" limit.
- b{3,10}_e{11,12,13,14}_f13: no STS.
- b{3,10}_e{13,14}_sts_f13: `use_rkg_sts = true`. Is a higher cap affordable with STS? Compare wall/rot
  with the non-STS twins.
- b{3,10}_e13_f16: the production floor.

bbot units: the code B is Heaviside-Lorentz (p_mag = B^2/2), so bbot_code = B_G/sqrt(4 pi). That gives
3 G = 0.846 and 10 G = 2.821. The field is a dipole with B_r = bbot cos(theta) (x1min/r)^3 at the inner wall.
cs_mhd_prod4's "bbot = 3 G" was code units, i.e. 10.6 G. Our 10 G arm matches it to within 6 %.
**`problem/bbot` in this package is in Heaviside-Lorentz CODE units, not Gauss: 3 G = 0.846, 10 G = 2.821**
(make_arm.sh converts BBOT_GAUSS for you). A physical-unit key `problem/bbot_gauss` (converted by the pgen,
startup print in both units, FATAL if both keys are set, bitwise with the old key) is on branch
`units-physical-inputs`, pending the user's merge decision; this package does not need it.

## Acceptance criterion for max_eta (state it in the results)
The acceptable cap is the LARGEST cap (closest to the physical eta) that meets all three conditions:
- (a) It is stable: no FATAL, NaN or dt collapse, and div B stays at round-off.
- (b) Its Ohmic dt costs at most ~2x the wall/rot of the eta ~ 0 arm with the same bbot and dfloor, or it
  is affordable with STS. Compare the STS twins.
- (c) Its energy-containing scales do not depend on the cap. Compare against the next-higher cap and the
  eta ~ 0 arm, at equal time (rot 301, 302):
  - the jet maximum and the equatorial u at 1e-3, 0.1 and 1 bar;
  - the day-night contrast at 0.1 bar;
  - Ohmic heating / L_abs;
  - t_drag = eta/v_A^2 against t_adv = pi r/|v_h| in the cells AT the cap. If t_drag(cap) >> t_adv
    there, raising the cap changes nothing dynamically.
  - The cap sweep in ana_eta.py measures, per pressure band, the volume fraction where the cap CREATES
    drag that the physical eta would not. This is t_drag(cap) < t_adv <= t_drag(uncapped). It should be
    about 0 for p > 1e-4 bar.
  - Only p > 1e-6 bar counts.
- Chaotic divergence: no noise member exists yet. Run one twin of one arm to measure the noise at equal time.
  For example, b3_e12_f13 at 2 ranks instead of 4: it is not bitwise, but it has the same physics.
  Differences below that twin spread are noise.

On the viper smoke state (3 G dipole, before any winding), the spurious-drag fraction for
p > 1e-4 bar was:
- 3-17 % at 1e11;
- at most 3 % at 1e12, and only at 1e-4 to 1e-2 bar;
- 0 at 1e13 and 1e14 (1e-6 to 1e-4 bar: 3.6 % at 1e13, 0.2 % at 1e14).

At 10 G, v_A^2 is 11x larger, so each row moves one decade: 1e13 at 10 G is about 1e12 at 3 G. The wound-up
field after 1-2 rotations decides.

## What to report (ana_eta.py per arm; paste the ana_eta.txt summary lines)
Per arm, report:
- stable? (FATAL/NaN/NOT-CONVERGED/dt COLLAPSE counts);
- dt mean/min and the limiter (Ohmic vs CFL, and where);
- wall per rotation and ms/cycle;
- f(eta at cap) per pressure band;
- the jet maximum and the day-night contrast;
- Q_ohm/L_abs;
- median t_drag/t_adv in the at-cap cells (capped / uncapped);
- the cap-sweep table;
- E_mag and eta_num.

eta_num is the estimated numerical diffusivity, (W_Lorentz - Q_ohm - dE_mag/dt)/int J^2 dV, from the
hst around each bin. It ignores the boundary Poynting flux, and the resolved J misses grid-scale
current, so it is an estimate. It needs the hst to cover +-0.1 rot around a bin.

Then give a recommended max_eta (and STS or not) for each bbot.
- Commit the RESULTS and the ana_eta.txt files (no bins or rst) on a new branch
  **bench-results-deltaai-w121mhd** and push it to fork.
- Do not merge it, and do not start any production.

## Known issues to watch
- **lhlld** (user rule 09-29 for dhj MHD):
  - Its low-Mach correction chi = min(1, max|u_n|/max c_f) has a MAGNETIC FLOOR. Where v_A dominates (the
    floored top, strong-field regions), chi -> 1 and it acts as plain hlld. The deep interior has beta >> 1,
    so there it is active.
  - `<mhd>/hlld_bx_zero_tol` is not set, so the historical default of 1e-4 is kept, as in prod4. This drops
    the rotational (**) states wherever beta > ~2e4, which is the whole deep interior here (beta ~ 1e9 at
    100 bar). That is stable at M_A >> 1 but more dissipative. The 1e-8 setting cured the Balsara-vortex
    instability at beta > 2e4 and M_A <~ 1 (lhlld_mhd.hpp note). It is not a scan variable.
  - Watch the low-beta top and the seam for checkerboards or dt collapse.
- The fresh start carries the hydro rot-300 state plus a pure dipole, so the first ~0.5 rot is field
  winding. Compare arms at equal time, not per cycle.
- The eta_eos cap sits in cold, metal-condensed gas: 0.01-1 bar, and on the night side up to about
  1e-4 bar. At 1e13, the fraction of cells at the cap per band was 20-27 % at 0.01-1 bar in the smoke state.

## USER DECISIONS 09-29 ~21:45
- **dfloor = 1e-13 for the MHD scan: YES** (the main scan as set up; keep the two 1e-16 check arms).
- **hlld_bx_zero_tol: keep the ORIGINAL value (1e-4, i.e. do not set the key).** No with/without pair needed.

## PACKAGE DOWNLOAD (09-29 ~22:10)
The package is on the fork, branch **data-w121-mhd-pkg** (data only, orphan history, 26.7 MB, no reference rst; every arm starts `-i` from remap.dat):
```bash
git fetch fork data-w121-mhd-pkg
git show fork/data-w121-mhd-pkg:w121_mhd_etamax_pkg_slim.tgz > w121_mhd_etamax_pkg_slim.tgz
md5sum w121_mhd_etamax_pkg_slim.tgz   # 7f02a2c1b8dd1be4830e2add7920d8fd
tar xzf w121_mhd_etamax_pkg_slim.tgz  # -> deltaai_pkg/
```
The full 128 MB tarball (with dhj.00600.rst) stays on viper only.
