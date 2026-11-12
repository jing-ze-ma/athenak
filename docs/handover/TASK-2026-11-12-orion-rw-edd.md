# TASK 2026-11-12 viper -> Orion: linear radiation-acoustic wave with the EDDINGTON closure (D = I/3) vs the analytic eigenmode (CPU)

## Why
The IDORT/VET grid (rc18_1110/RW_AMP1E6_1111.md, TASK-2026-11-12-orion-rw-p333) scores the code against vetqs_d, a
reference whose closure model may not be exactly what the code solves (P 3.33 tau 10 converges to an omega_inf off
vetqs_d by 6.5e-4 rel.). Here the code's closure is fixed to the Eddington tensor, closure = eddington (rad_m1.cpp:
chi = 1/3 exactly, D off-diagonal 0; no rays, so implicit_flux_beam = halfrange leaves every face central, rad_m1_implicit.cpp
"without ray data ... the faces stay central unless implicit_hr_model"), i.e. the two-moment Eddington system of Jiang,
Stone & Davis 2012 (ApJS 199, 14, Sect. VI.1.1 / App. B). The reference is then exact: no model mismatch. Expect clean
2nd order in every cell incl. P 3.33 tau 10; any cell that is not = a scheme defect. No code change needed.

## Builds (reuse the worktrees of TASK-2026-11-12-orion-rw-amp1e6; same cmake line there: Release, MPI on, build_none)
- PROD arm: 7d233ac30 (fork/vgdhr-1110) -> <wt_rwa6_prod>/build_none/src/athena
- FINAL arm: 84cd20888 (fork/final-1111; rad_m1 identical to fork/final2-1111 156d7ccd0) -> <wt_rwa6_fin>/build_none/src/athena

## Cases (docs/handover/scripts/rw_edd; runner = ../rw_amp1e6, keep the two directories siblings)
18 gated cells: P_rad/P_gas 0.00333 (tau 0.1, 10, 1000), 0.0333 and 0.333 (tau 0.01, 0.1, 1, 10, 100, 1000), 3.33 (tau 0.1,
10, 1000) x N 16, 32, 64, 128, 256, nt = 2N, x1, transport implicit, hesdirk2, vimp, ref = edd; 90 cases per arm.
Every case: amp 1e-6, implicit_tol 1e-12, implicit_lin_tol 1e-13, implicit_lin_cnorm 1e-13, implicit_tol_rel 0,
dt_retry false (drive.py forces it). PROD = arm PROD19 keys; FINAL = arm TIME keys + implicit_physflux_zold = fc +
implicit_hr_cmem_rst = true. c = 10, chat/c = 1, gamma 5/3, kappa_P = kappa_E = kappa_F, no scattering.
```
cd <run area>/rw_edd           # copy of docs/handover/scripts/rw_edd, with ../rw_amp1e6 next to it
export OMP_NUM_THREADS=1 RW_PKF=$PWD/../rw_amp1e6/prodkeys_all.txt
# check the input first: closure = eddington, dt_retry = false, implicit_tol = 1.0e-12, radwave_amp = 1e-06
# smoke (N16, 18 cases per arm; no FATAL expected -- if one FATALs on a PROD19/TIME key with eddington, report it):
RW_EXE=<wt_rwa6_prod>/build_none/src/athena RW_RUNS=$PWD/runs_PROD  python3 ../rw_amp1e6/drive.py run smoke_PROD.txt  --np 18 > smoke_PROD.out 2>&1
RW_EXE=<wt_rwa6_fin>/build_none/src/athena  RW_RUNS=$PWD/runs_FINAL python3 ../rw_amp1e6/drive.py run smoke_FINAL.txt --np 18 > smoke_FINAL.out 2>&1
# full (the N16 cases are cached from the smoke; lists are largest-N first):
RW_EXE=<wt_rwa6_prod>/build_none/src/athena RW_RUNS=$PWD/runs_PROD  python3 ../rw_amp1e6/drive.py run cases_PROD.txt  --np 18 > cases_PROD.out 2>&1 &
RW_EXE=<wt_rwa6_fin>/build_none/src/athena  RW_RUNS=$PWD/runs_FINAL python3 ../rw_amp1e6/drive.py run cases_FINAL.txt --np 18 > cases_FINAL.out 2>&1 &
wait
```
Cost: no VET sweeps with eddington, so cheaper than the vetf grid (vetf N256 nt512 was ~23 s on a viper core; the thin
tau 0.01 vetf cost was the sweep passes, absent here). A FATAL is a result: report it, do not re-enable dt_retry.

## Analytic reference (verified on viper)
`python3 ref_edd.py`: the acoustic root of radwave_disp.matrix (the lab-frame linearised system as the code writes it:
gas momentum + kappa (F - 4/3 E0 v)/c, gas energy (gamma-1) c kappa (E - 4 a T0^3 T), radiation E/F with the c rho kappa
exchange and drag, P1 = E1/3). Limits check exactly: tau -> 0 c_ph = sqrt(gamma T0) = 1.29099445; tau -> inf c_ph =
c_mix / sqrt(1 + 4 E0 / (3 rho c^2)) (P 3.33: 2.29137806). The lab form matters at large P: eddc (comoving F0, no
-(4/3) E0 dv/dt) differs by 8.4 % at P 3.33 tau 10 (edd 13.70132 - 4.08637 i, eddc 14.63663 - 4.83298 i); earlier VET
runs matched the lab form (vetqs_d) to 6e-5, so ref = edd. analyze_edd.py also prints |d omega| vs eddc to confirm.

## Return (NOTE-2026-11-1x-orion-viper-rw-edd.md on this branch; NOTE only, no run trees)
1. `RWE_ROOT=$PWD python3 analyze_edd.py` (orders p(rho), p(E) for each N doubling, mean-removed p(E) 128->256,
   N256 L1, |omega - omega_edd|/|omega| and vs eddc) -- paste whole.
2. `grep -c "ok w" cases_*.out; grep -h -E "FAILED|FATAL" cases_*.out` and the "ok w=... ref=..." line of every N256 case.
3. `python3 ref_edd.py` output; wall time per arm; hostname / cores used.
