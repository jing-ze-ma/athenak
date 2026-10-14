# NOTE Caltech: bsgli QUEUED (TASK-2026-10-14-caltech-bsgli)

**No smoke/gate: user decision 10-14, production submitted directly.** run.cfg BLC_GATE=0 (bsgli_hook.sh line 21 checks
VERIFIED only when BLC_GATE=1). No other change from the TASK.

- Binary /resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_8e018272_nofma, md5 4a92728d016609effe6d4586fb2e55bf,
  git rev-parse HEAD 8e0182723da23d8c864c8b079bd24e75689f2d28 (= origin/fin4-linhao-bc-8e018272); -ffp-contract=off, HOPPER90, MPI.
- Deck bsgli_1014/files/bsgLI_bc_caltech.athinput (read-only), md5 d77f9a102c2fcc4db73f97422cee7d56 (template 532fc7af, @FILES@ sed only).
- linhao_profile_full.txt md5 41edc225; ic_ma2026_f0 246a666e, rosseland 450bc0c1, planck e18ca805; all MD5SUMS pass.
- Layout 1 node x 2 H200, 2 ranks, LAUNCH = e1_1013/e1launch.sh, --exclude=hpc-sm-01-09.
- Run dir /resnick/groups/carnegie_poc/jingze/bsgli_1014/prod (fresh). TLIM 4.96e6, STOP_NC=1.
- Links (-J bsgLIc, 12 h): 4487415 (no dependency), 4487610, 4487611, 4487612, 4487613 (afterany chain).
- Start estimate: sbatch --test-only 2026-10-13T23:49 (hpc-sm-01-05); squeue estimate for 4487415 2026-10-14T11:40.

Next: STARTED when link 1 runs, then per-link reports as in TASK section 6.

## STARTED
- Link 1 4487415 RUNNING on hpc-sm-01-05 since 2026-10-13T21:02:01, fresh t = 0 (no rst in prod/rst).
- he_ic_linhao line: linhao_profile_full.txt (327682 rows, r 1.14819e+12 .. rmax 3.65147e+12), T(rmax) 63364.2 K, rho(rmax) 2.64193e-11.
- he_bc_linhao ghost lines printed (ghost 0: rho 6.462e-05, T 1.0963e6, F_r 2.330e13).
- vet_gd_halo_ipc: ON (node-local CUDA IPC). dt 87.5 s from cycle 50; ~0.86 s/cycle over cycles 50-150.
- NOT ADMISSIBLE (redone with BE) at cycles 3 and 4 only so far; 0 FATAL / DIVERGED / NON-CONVERGED.

## EARLY NUMBERS (link 1 4487415, t 8.1e4 s, cycle ~970; USER asked to send) -- question: do Raven bsgLI / Viper bsgLIv show the same?
Clean: 0 FATAL / NaN / DIVERGED / NON-CONVERGED; NOT ADMISSIBLE only 4 (cycles 3-4, start); Picard 3-4; dt 87.5 -> 83.33 s at ~4.4e4 s;
~0.9 s/cycle (wall 794 s for 900 cycles incl. init). 8 bins + rst 00000.

user.hst (wedge values):
| t [s] | dt | L_bot | L_mid | L_top | L_in | KE_int | Etot | Picard |
|---|---|---|---|---|---|---|---|---|
| 0 | 87.50 | 4.651e37 | 4.651e37 | 0 | 4.651e37 | 0 | 1.56967e46 | 0 |
| 1050 | 87.50 | 4.586e37 | 4.051e37 | 1.004e38 | 8.546e38 | 1.32e40 | 1.56912e46 | 8 |
| 5075 | 87.50 | 4.860e37 | 4.428e37 | 5.004e37 | 2.738e39 | 5.32e40 | 1.56928e46 | 3 |
| 14000 | 87.50 | 5.541e37 | 4.410e37 | 4.578e37 | 3.735e39 | 1.36e41 | 1.57003e46 | 3 |
| 29050 | 87.50 | 6.735e37 | 4.470e37 | 4.591e37 | 3.517e39 | 1.67e41 | 1.57152e46 | 3 |
| 44000 | 83.33 | 7.548e37 | 4.536e37 | 4.616e37 | 2.732e39 | 1.37e41 | 1.57276e46 | 3 |
| 59000 | 83.33 | 7.978e37 | 4.531e37 | 4.612e37 | 2.026e39 | 1.20e41 | 1.57364e46 | 3.4 |
| 81000 | 83.33 | 8.150e37 | 4.411e37 | 4.505e37 | 1.374e39 | 1.24e41 | 1.57441e46 | 3.8 |

Observations (not judged; may be the intended he_bc_linhao behaviour):
- L_mid / L_top stay at the IC L (4.4-4.6e37) after the start transient.
- L_bot rises 4.65e37 -> 8.15e37 (x1.75), flattening.
- L_in (the 8e018272 diagnostic column) peaks 3.8e39 at ~1.7e4 s (~80x L_bot), now 1.37e39 and falling.
- Etot rises ~7e38 erg/s early, ~5.6e38 erg/s over 2.9e4-8.1e4 s (~12x L_top): energy enters through the inner boundary
  (hydro ghosts / ghost bath) faster than radiation leaves at the top.
- KE peaks 1.7e41 at ~2.5e4 s (the 5 % rho perturbation), now ~1.2e41.
