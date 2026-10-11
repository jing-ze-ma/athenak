# NOTE gdface (gdface-1011): face-consistent inner diffusion intensity, now the default everywhere

Branch `deepatm-1011` in /raven/u/jinma/ATHENAK/wt-deepatm. Local commits only (not pushed):
- **3b50c192**: vet_sc `vet_bc_diff = face` becomes the DEFAULT. Opt-in face start for vet_gd, vet_col_lat and the vet_col core.
- **e84d2513**: `face` becomes the DEFAULT under closure = vet_col too (vet_gd, vet_col_lat, vet_col core).

Work dir: /raven/ptmp/jinma/deepatm_1011/gdface.

## 1. What changed

There is one key, `<rad_m1>/vet_bc_diff = face | cell`.
- The closures are mutually exclusive, so a single key is unambiguous. It is read in VetInit for vet_sc and in the vet_col branch of the constructor for vet_col / vet_gd / vet_col_lat.
- One knob reproduces any old run.
- The default is `face`. A **restart whose file lacks the key keeps `cell`**, the same convention as vet_bc_bath and vet_col_surface_q. The resolved value is echoed in the closure banner.

| path | cell (old) | face (new default) |
|---|---|---|
| vet_sc bottom (all 4 sweep paths) | eps = E_0 + 3F.Omega/c | E_f = E_0 + (S_f - S_0) (deepatm f0b6b8a4) |
| vet_gd inner face of shell icut (rad_m1_vetgd.cpp ~2000) | S_a + mu(-dS/dtau) | (1.5 S_a - 0.5 S_b) + mu(-dS/dtau) |
| vet_col_lat VL_START (rad_m1_vetlat.cpp ~713) | same as vet_gd | same as vet_gd |
| vet_col core rays without vet_col_order2 | E_0 + 3 F_0 mu/c | E_0 + 0.5(S_0 - S_1) + 3 F_0 mu/c (as vet_sc) |
| vet_col core rays with vet_col_order2 (default on sp) | already face (1.5E_0 - 0.5E_1, face flux) | unchanged |

In every case the half-cell carry already started from the face source S_f = 1.5 S_0 - 0.5 S_1. The new start makes the boundary intensity refer to the same point.

GPU determinism and halo:
- The change is pure per-thread arithmetic, with no new data, reductions or messages. vet_gd_twin_det and the compact halo are untouched.
- face is run-to-run bitwise on 4 A100: 18 bin files + 2 hst identical.
- The halo_compact = 1 vs 2 check (hc1 arm) runs in the combined job, see section 4.
- vet_gd_halo_compact = 0 is FATAL with vet_gd_twin_fuse + MPI. That is a pre-existing constraint, not a regression.

## 2. Bitwise status

| check | result |
|---|---|
| CPU atm ladder 32..2048: new default tJe / cen vs old f0b6b8a4 tJeD / cenD | 28/28 runs, 0 differing data |
| CPU atm: new tJeC / cenC (vet_bc_diff = cell) vs base 42aa6bbc tJe / cen | bitwise (in the 28) |
| CPU 1-D ladders: pulse_k0.128 + rw_t10, be + hesdirk2, X/T/C, 60 runs vs old tJeD | 0 differing |
| G1 (cfl 1e2, 1e4), G3 (c1, c10), old G5 (4 kappa): tJe vs old tJeD; tJeC vs base tJe | 16/16 bitwise |
| beams ba0, xb20 (tJe vs old tJeD) | bitwise (Picard 4.19 / 4.86, 0 non-converged) |
| vet_col atm 128/512: e84d default == 3b50 face; e84 cell == 3b50 cell | bitwise |
| GPU AG Car A 50 cycles: e84d default == 3b50 explicit face; e84 cell == 3b50 cell (== the old code) | bitwise (18 bin + 2 hst) |
| GPU AG Car A face, repeated | bitwise (twin_det) |

The e84d2513 binaries are bin/athena_{cpu,he_a100}_wipB, built from src_wipB, which is byte-identical to the commit's rad_m1 sources.

Numbers:
- atm orders, tau > 10, 32..2048:
  - tJe (face): 1.99 2.00 2.00 2.00 2.00.
  - tJeC (cell): 2.01 2.04 1.77 1.49 1.47.
  - cen: 2.04 2.08 2.20 2.36 2.18. cenC: 2.06 2.13 1.91 1.57 1.38.
- G1 Hopf L1: 6.994e-4 (face) vs 7.080e-4 (cell).
- G3: L1 0.0139 / 0.0142 PASS.
- Fixed G5 (g5fix pulse2d_g5, run in g5/), all 12 arms PASS, the same ratios as the g5fix baseline:
  - cen: 0.999609 0.999977 1.000000 1.000000.
  - tJe_be: 0.999581 0.999972 1.000000 1.000000.
  - tJe_h2: 0.999565 0.999874 1.000000 1.000000.
- vet_col (plane-parallel atm, nmu 4, arms vcolC / vcolD, run in vcol_run/):
  - It is second order in every band both ways. The orders are identical to the printed digits.
  - The face-cell difference is <= 1e-11 relative, i.e. at the solver tolerance.
  - So vet_col had the mismatch, and the fix is applied for consistency, but it is neutral: vet_col's uniaxial f_K is insensitive to a half-range deficit at depth (the D shift is proportional to (D - 1/3)), and vet_col's H is not used.

**To reproduce any run made before e84d2513 / 3b50c192 bitwise, add `vet_bc_diff = cell` to `<rad_m1>`**. A restart from an old restart file already resolves to cell.

## 3. Production impact (Q3)

GPU A/B on 1 node x 4 A100 (job 31065032), cell vs face, from t = 0:
- AG Car A: agcar_A_xtf, 50 cycles.
- He giant (vgdfuse heg) and BSG (vgdfuse bsg reduced to 4 blocks): 30 cycles each.
- Dumps at the last cycle. Analysis: ana_ab.py, RESULTS_ab.txt.

Where the boundary sits:
- vet_gd starts at shell icut, where tau from the top reaches vet_col_lat_taucut = 30. Below icut the tensor is isotropic (vet_gd_replace).
- The domain bottom is optically enormous, so its own dtau is irrelevant: no formal solution starts there.
- The relevant cell is the icut shell. Its dtau_cell is 4-10, so the face-cell offset 0.5(S_a - S_b) is O(1) of 3H dtau. That is large for H right there, but small relative to J ~ 3H tau (tau ~ 35-42).

| model | innermost cell dtau (domain base, IC + kappa_R) | icut shell / r | dtau_cell(icut) | tau at its inner face | photosphere shell (tau 2/3) |
|---|---|---|---|---|---|
| AG Car A | 4.9e4 (tau_base 3.8e6) | 248 / 379.4 Rsun | 4.2 | 35.6 | 268 (388.3 Rsun) |
| He giant | 6.9e6 (tau_base 2.0e8) | 264 / 61.6 Rsun | 10.5 | 42.4 | 271 (61.9 Rsun) |
| BSG | 6.1e5 (tau_base 5.8e6) | 151 / 51.2 Rsun | 7.8 | 41.1 | 163 (53.1 Rsun) |

Bottom-shell flux ratio of the formal solution, H_gd/(F/c), shells icut..icut+3 (exact ~1 in the deep layer):

| model | cell | face |
|---|---|---|
| AG Car A | 0.903 1.004 1.000 1.000 | **0.976** 1.000 1.000 1.000 |
| He giant | 0.349 0.732 0.978 1.121 | 0.351 0.732 0.978 1.121 |
| BSG | 0.957 1.001 1.000 1.000 | **0.973** 1.001 1.000 1.000 |

He giant:
- Its first shells are wrong for a different reason, and that error dominates: dtau_cell 10.5 with the linear-in-tau source over steep S.
- The face start changes that only by 0.002, because the half-cell carry attenuates the boundary intensity by exp(-5).
- Fixing it needs resolution at icut or a larger taucut, not the boundary intensity.
- H_gd is not used by these production inputs (no thin switch). It enters only through D.

D_rr at icut (exact 1/3 in the deep layer):

| model | cell | face |
|---|---|---|
| AG Car | 0.333006 | 0.333308 |
| BSG | 0.333230 | 0.333305 |
| He giant | 0.333321 | 0.333323 |

Field differences, face vs cell. The metric is the per-shell max over (theta, phi) of |a-b| / max|b|.

| | AG Car A (50 cyc) | He giant (30) | BSG (30) |
|---|---|---|---|
| E at icut (largest anywhere) | 8.9e-4 | 9.1e-6 | 2.2e-4 |
| D_rr at icut | 9.2e-4 | 6.6e-6 | 2.2e-4 |
| F1, any shell | <= 8.4e-7 | 2.5e-4 (icut+1) | <= 5.4e-7 |
| dens at icut | 2.4e-4 | 2.6e-7 | 2.0e-5 |
| dens at the photosphere | 8.6e-8 | 0 | 1.1e-7 |
| v_r at the photosphere / global max\|v\| | 1.2e-7 | 4.0e-8 | 1.2e-6 |
| v_r max / global max\|v\| | 6.7e-3 (r 407 Rsun, rho 1e-16 floor envelope) | 5.1e-6 | 1.3e-4 (icut-1) |
| hst L_top, max over cycles | 4.0e-6 | 1.2e-8 | 3.6e-6 |
| hst Mdot_top | 0 both (no outflow yet) | 0 both | 0 both |
| s/cycle cell / face | 0.412 / 0.400 | 1.177 / 1.180 | 0.330 / 0.324 |
| Picard mean (max), cell / face | 3.9 (9) / 3.9 (9) | 9.5 (12) / 9.5 (12) | 2.70 (6) / 2.67 (6), 0 non-converged |

The AG Car v_r maximum sits in the near-floor envelope above the photosphere (rho ~ 1e-16 g/cc, 407 Rsun):
- That gas has almost no inertia and responds to any perturbation of the radiation field.
- In the dense layers the difference stays at 1e-7 or below.
- The noise-reference arm tol9 (implicit_tol 1e-9 instead of 1e-8, same binary) runs in the combined job, see section 4.

### Plain answer to Q3

**No, not for the dynamics.**
- In all three production models the formal solution starts at tau ~ 35-42, only 7-20 cells below the photosphere.
- The face fix changes only the first one or two shells there:
  - H_gd moves toward the exact value (AG Car 0.90 -> 0.98, BSG 0.96 -> 0.97, He giant unchanged at 0.35).
  - D_rr moves toward 1/3 by ~3e-4.
  - E changes locally by <= 9e-4.
- Everywhere else the change is <= 1e-6, including the photosphere, L_top (<= 4e-6) and Mdot_top (0).
- Cost and Picard are unchanged.
- The effect is small because D = K/J is insensitive to a half-range deficit at depth (the shift is proportional to D - 1/3), and these inputs do not use H_gd.
- It is not zero, because dtau_cell(icut) = 4-10 makes the offset O(1) for H in the first shell.

The large deep-atm effect of deepatm-1011 (order 1.5 -> 2) is specific to vet_sc, whose H_sc enters the solution, and to resolved deep layers.

Production runs should switch to the default at their next fresh start. A continuing restart keeps `cell` automatically and stays bitwise.

## 4. Pending (combined GPU job of the icrad agent, see COMBINED_JOB.txt)

- **tol9**: agc cell with implicit_tol 1e-9. The noise reference for the A/B differences.
- **hc1**: agc face with vet_gd_halo_compact = 1. It should be bitwise to face (halo_compact 2).

At 04:47 the job was still pending (Priority), with an estimated start of 10:40.

To finish when it is done (S = runs/j.31066055, R = runs/j.31065032):
- `bash ../bitcmp.sh $R/agc_face $S/agc_face_hc1`: halo exactness. Expect 0 differing.
- `python ana_ab.py $R/agc_cell $S/agc_cell_tol9 tol9`: the size of a solver-tolerance perturbation, to compare with the cell-vs-face table above. Mainly the v_r in the floor envelope at r ~ 407 Rsun.

## Files

- Scripts: build.sh, cpu_val.sh, gpu_ab.sh (`sbatch gpu_ab.sh <bin> "<case:arm[:k=v,..][:tag][:bin]> ..." NA NB`), ana_ab.py, tau_est.py, mkinp.py.
- Inputs: inp/{agc,heg,bsg}.athinput. These are the production inputs with the dump set hydro_w, m1, m1_vet, m1_fs, m1_face and m1_thsw at the last cycle.
- Runs:
  - cpu_run/ (atm, ladders, gates, beams).
  - vcol_run/ (vet_col atm).
  - g5/run.
  - vcolchk/, smoke*/ (checks).
  - runs/j.31065032 (GPU A/B).
- Results: RESULTS_ab.txt, g5/g5_results.txt.
- Binaries:
  - bin/athena_{cpu,he_a100}_3b50c192 (md5 7b554e06 / 63785815).
  - bin/athena_{cpu,he_a100}_wipB (= e84d2513, md5 d26a2d82 / b753c4a0).
- Changes to the deepatm cpu bundle: battery_raven.sh A[tJeC], A[cenC]; od.py ARMS tJeC, vcol, vcolD, vcolC. Backups are *.bak_gdface.
- Known diagnostic quirk (not changed): with closure = vet_col, an `m1_vsc` tab output block suppresses the final time-based dump. od.py's atm template has one, so vcol_run strips it.
