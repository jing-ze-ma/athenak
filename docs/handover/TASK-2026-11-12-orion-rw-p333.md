# TASK 2026-11-12 viper -> Orion: IDORT linear wave, cell P_rad/P_gas 3.33 tau_lambda 10 -- why the order is lost (CPU)

## Why
In the amp 1e-6 grid (viper rc18_1110/RW_AMP1E6_1111.md) the cell P 3.33, tau 10 is the one non-thin cell that loses
order: p(rho) 1.12, p(E) 1.37 at N 128 -> 256, PROD and FINAL alike, amp 1e-5 and 1e-6 alike. Viper's existing runs
(FINAL binary 84cd20888, all with dt_retry off) already show:

| run (nt) | w_re | w_im |
|---|---|---|
| reference vetqs_d (nmu 4 x nphi 8) | 13.346366 | -3.660631 |
| FINAL amp 1e-6, nt = 2N: N64 / N128 / N256 | 13.325986 / 13.340693 / 13.344398 | -3.663332 / -3.663039 / -3.663014 |
| amp 1e-5 tol 1e-11, N256 time ladder nt 512 / 1024 / 2048 / 4096 | 13.344398 / 13.344047 / 13.343960 / 13.343939 | -3.663014 / -3.662965 / -3.662953 / -3.662951 |
| amp 1e-5, fixed c dt/dx 0.588: N32 / 64 / 128 / 256 | 13.244826 / 13.319473 / 13.338975 / 13.343960 | -3.660900 / -3.662290 / -3.662794 / -3.662953 |

i.e. dt converges at 2nd order (differences 4.9e-5, 1.2e-5, 2.4e-6) and dx at ~2nd order, both to
omega_inf ~ (13.3456, -3.6630) != reference (13.34637, -3.66063): the damping rate is off by ~ -2.4e-3 (6.5e-4
relative, = the -6e-4/period amplitude floor). nmu 8 / nphi 16 changes w_im by 1e-5 only. Hypothesis: the code
converges to a slightly different continuum model than the vetqs_d reference. The runs below (amp 1e-6, tight tol,
dt_retry off -- user rule 11-12) confirm the limits and find which ingredient moves omega_inf.

## Build
FINAL arm only: commit 84cd20888 (fork/final-1111; rad_m1 identical to fork/final2-1111 156d7ccd0). The worktree
wt_rwa6_fin of TASK-2026-11-12-orion-rw-amp1e6 can be reused (build_none, MPI on, Release).

## Runner (docs/handover/scripts/rw_p333 + rw_amp1e6 next to it; copy BOTH directories to the run area, keep them siblings)
```
cd <run area>/rw_p333
export OMP_NUM_THREADS=1 RW_PKF=$PWD/../rw_amp1e6/prodkeys_all.txt RW_EXE=<wt_rwa6_fin>/build_none/src/athena RW_RUNS=$PWD/runs
python3 ../rw_amp1e6/drive.py run cases_base.txt --np 6  > cases_base.out 2>&1
python3 ../rw_amp1e6/drive.py run cases_a.txt    --np 29 > cases_a.out 2>&1
python3 drive_diag.py          run cases_diag.txt --np 2  > cases_diag.out 2>&1
```
(the three can run concurrently; 37 cores). Check first: `python3 drive_diag.py input cases_diag.txt | grep -E
"dt_retry =|^implicit_tol|lin_cnorm|radwave_amp|variable"` -> dt_retry = false, implicit_tol = 1.0e-12, tol_rel 0,
cnorm 1.0e-13, amp 1e-06, outputs m1_vet and m1_vsc.

| list | cases | purpose |
|---|---|---|
| cases_base.txt | nt = 2N, N 16..256 and N1024 nt2048 | same-build baseline + the N->inf limit |
| cases_a.txt 1-4 | N512 nt1024; N256 nt 1024/2048/4096 | (a) time ladder at N256 |
| cases_a.txt 5-8 | nt 4096 fixed, N 32/64/128/512 (N256 nt4096 is line 4) | (b) self-convergence in space at fixed small dt |
| cases_a.txt 9-29 | N 128/256/512 (nt 2N) x 7 arms: implicit_flux=central, vet_source=gas, time2_vet_sc=rebuild, implicit_closure_lag=pass, vet_scatter=false, implicit_enthalpy=upwind, implicit_physflux=false | (d) closure update within the step (rebuild, pass) and which ingredient moves omega_inf |
| cases_diag.txt | N128, N256 at nt 4096 with m1_vet / m1_vsc dumps | (c) the code's VET perturbation vs the vetqs_d closure |

Cost on one viper core: N256 nt512 23 s; the longest, N512 nt4096 and N1024 nt2048, ~6-10 min. A FATAL is a result
(report it; e.g. physflux=false may refuse a FINAL key), do not re-enable dt_retry.

## Return (NOTE-2026-11-1x-orion-viper-rw-p333.md on this branch)
1. `python3 reana.py cases_base.txt runs; python3 reana.py cases_a.txt runs` (N nt c dt/dx w_re w_im L1rho L1E t_end
   arm) -- paste both outputs whole.
2. `python3 diag_ana.py runs/<N256 nt4096 diagvet dir>` and the same for N128 (per dump: R = E0 dD11/(qT T1), the code's
   E1/T1 and T1/rho1 relative to the vetqs_d eigenvector) -- paste whole.
3. `python3 variants.py` output; wall per list. Commit the NOTE only (no run trees).
