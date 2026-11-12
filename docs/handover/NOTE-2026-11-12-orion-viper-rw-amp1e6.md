# NOTE 2026-11-12 Orion -> viper: rw-amp1e6 (N256 tolerance check + addendum grid cases)

Builds: PROD 7d233ac30 (athena md5 423209affc92f7e157c24c3194a97357), FINAL 84cd20888 (8ec46383d524802271739678056f13c4);
Release, MPI, -ffp-contract=off. One rank per case, OMP 1, Orion p.shared (orion001 / orion003). The pre-check passed for every list:
dt_retry = false, implicit_tol 1e-13 (tt) / 1e-12 (grid), lin_tol and cnorm 1e-14 / 1e-13, tol_rel 0, radwave_amp 1e-06.
analyze.py is the addendum version; n128_ref/runs_* were copied in (cp -n) before the addendum runs.

## Verdict
- **All 8 tolerance-check cases (orion_tt256_*) and all P 0.0333 / 0.333 addendum cases end "ok"**: Picard mean 3 (max 3), NON-CONVERGED 0 in every case.
- **Tolerance-converged:** tt (1e-13/1e-14) vs grid (1e-12/1e-13) at N256 on the same Orion build:
  - P 0.0333 tau 0.01: rho, E and E-mrm equal to 4 digits (PROD and FINAL);
  - P 0.333 tau 0.01: equal for PROD; FINAL E 2.745e-4 vs 2.749e-4 (0.15 %).
  So far inside 1 %. The tt rows for P 0.0333 tau 0.1 and P 3.33 tau 10 are for your grid comparison.
- **Thin tau 0.01 orders 128 -> 256** (now complete with n128_ref): p(rho) 2.00 / 2.07; p(E) raw 1.91, mean-removed 1.93-1.94; FINAL/PROD L1 at N256 1.000 / 0.999.
- **P 33.3 tau 0.01 N256 (info only).**
  - PROD FATALs: hesdirk2 stage 2 not admissible at cycle 46, then a non-finite hydro state at cycle 47.
  - FINAL does not FATAL, but the solution diverges: w = 7.566 + 39.25i (growing; ref 36.25 - 0.082i), L1 rho 3.2e5, Picard mean 18.9 / max 47,
    682 solves instead of 1023. **Treat the FINAL "ok" as a failed (unstable) run, not a result.**

## RWA_ROOT=$PWD python3 analyze.py
```
| P | tau | arm | p rho | p E raw | p E mean-rm | p E raw 64->128 | L1 rho 256 | L1 E 256 | L1 E-mrm 256 | mean off 256 (A|dE|) | done N |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.0333 | 0.01 | PROD | 2.00 | 1.91 | 1.94 | 1.89 | 6.939e-04 | 2.799e-04 | 2.745e-04 | 7.920e-05 | 16,32,64,128,256 |
| 0.0333 | 0.01 | FINAL | 2.00 | 1.91 | 1.93 | 1.89 | 6.939e-04 | 2.801e-04 | 2.748e-04 | 7.914e-05 | 16,32,64,128,256 |
| 0.0333 | 0.01 | PROD_tt | 2.00 | 1.91 | 1.94 | 1.89 | 6.939e-04 | 2.799e-04 | 2.745e-04 | 7.920e-05 | 64,128,256 |
| 0.0333 | 0.01 | FINAL_tt | 2.00 | 1.91 | 1.93 | 1.89 | 6.939e-04 | 2.801e-04 | 2.748e-04 | 7.914e-05 | 64,128,256 |
| 0.0333 | 0.1 | PROD_tt | - | - | - | - | 6.005e-04 | 2.917e-04 | 2.908e-04 | 3.455e-05 | 256 |
| 0.0333 | 0.1 | FINAL_tt | - | - | - | - | 6.005e-04 | 2.917e-04 | 2.908e-04 | 3.455e-05 | 256 |
| 0.3333 | 0.01 | PROD | 2.07 | 1.91 | 1.93 | 1.89 | 5.956e-04 | 2.747e-04 | 2.708e-04 | 6.276e-05 | 16,32,64,128,256 |
| 0.3333 | 0.01 | FINAL | 2.07 | 1.91 | 1.93 | 1.89 | 5.956e-04 | 2.745e-04 | 2.706e-04 | 6.302e-05 | 16,32,64,128,256 |
| 0.3333 | 0.01 | PROD_tt | - | - | - | - | 5.956e-04 | 2.747e-04 | 2.708e-04 | 6.276e-05 | 256 |
| 0.3333 | 0.01 | FINAL_tt | - | - | - | - | 5.956e-04 | 2.749e-04 | 2.710e-04 | 6.276e-05 | 256 |
| 3.333 | 10 | PROD_tt | - | - | - | - | 8.529e-05 | 1.075e-04 | 1.075e-04 | 9.356e-08 | 256 |
| 3.333 | 10 | FINAL_tt | - | - | - | - | 8.513e-05 | 1.073e-04 | 1.073e-04 | 9.356e-08 | 256 |
| 33.33 | 0.01 | FINAL | - | - | - | - | 3.158e+05 | 1.413e+03 | 1.382e+03 | -1.367e+02 | 64,256 |

FINAL/PROD L1 at N256 (rho, E):
  P 0.0333 tau 0.01: 1.000 1.000
  P 0.3333 tau 0.01: 1.000 0.999
```

## Per case: Picard / NC counters (log tail) and cpu time (= wall, 1 rank)
```
runs_ttPROD | P0.0333333_t0.01_N256_nt512 PROD19 | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 5.430652e+03 | ok
runs_ttPROD | P0.0333333_t0.1_N256_nt512 PROD19 | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 5.690673e+02 | ok
runs_ttPROD | P0.333333_t0.01_N256_nt512 PROD19 | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 5.542942e+03 | ok
runs_ttPROD | P3.33333_t10_N256_nt512 PROD19 | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 2.388452e+01 | ok
runs_ttFINAL | P0.0333333_t0.01_N128_nt256 TIME | solves 5.110000e+02 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 1.204349e+03 | ok
runs_ttFINAL | P0.0333333_t0.01_N256_nt512 TIME | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 5.452958e+03 | ok
runs_ttFINAL | P0.0333333_t0.1_N256_nt512 TIME | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 5.716494e+02 | ok
runs_ttFINAL | P0.333333_t0.01_N256_nt512 TIME | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 5.545518e+03 | ok
runs_ttFINAL | P3.33333_t10_N256_nt512 TIME | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 2.496833e+01 | ok
runs_PROD | P0.0333333_t0.01_N256_nt512 PROD19 | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 4.887217e+03 | ok
runs_PROD | P0.333333_t0.01_N256_nt512 PROD19 | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 4.992740e+03 | ok
runs_PROD | P33.3333_t0.01_N256_nt512 PROD19 |  | cpu_s  | FATAL
runs_FINAL | P0.0333333_t0.01_N128_nt256 TIME | solves 5.110000e+02 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 1.206708e+03 | ok
runs_FINAL | P0.0333333_t0.01_N256_nt512 TIME | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 4.901698e+03 | ok
runs_FINAL | P0.333333_t0.01_N128_nt256 TIME | solves 5.130000e+02 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 1.233574e+03 | ok
runs_FINAL | P0.333333_t0.01_N256_nt512 TIME | solves 1.023000e+03 picard mean 3.000000e+00 max 3.000000e+00 NC 0.000000e+00 | cpu_s 5.001300e+03 | ok
runs_FINAL | P33.3333_t0.01_N256_nt512 TIME | solves 6.820000e+02 picard mean 1.892962e+01 max 4.700000e+01 NC 0.000000e+00 | cpu_s 5.398502e+03 | ok
```

## P 33.3 PROD FATAL (log tail)
```
<rad_m1> hesdirk2 stage 2 NOT ADMISSIBLE at cycle 46 (converged=1 min(E,T)=-1.000000e+00): the step is redone with backward Euler
### non-finite hydro state at cycle=47 time=1.581153e-02 cell (m,k,j,i) = (0,0,5,139) gid = 0 rank 0 r=5.371094e-01 rho=nan T=-nan p=-nan v=(nan,nan,nan)
### FATAL ERROR in /orion/ptmp/jinma/rwamp6_1112/src_7d233ac30/src/src/mesh/mesh.cpp at line 1547
dt = -8.000000e-01 at cycle 47 time 1.581153e-02: non-finite or non-positive time step (rank 0)
```

## results.json
Every Orion-run case is at docs/handover/scripts/rw_amp1e6/results_amp1e6_orion/<runs set>/<case dir>/results.json: tt PROD 4, tt FINAL 5 (incl. tt128), grid PROD 2, grid FINAL 5 (incl. N128 x2 and P 33.3).
The P 33.3 PROD FATAL has none.
