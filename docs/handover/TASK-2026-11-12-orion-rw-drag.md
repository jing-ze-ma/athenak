# TASK 2026-11-12 viper -> Orion: rw-drag, localise the ~3e-4 drag deficit at P 3.33 tau 10 (CPU)

## Why (viper rc18_1110/RW_P333_TERM_1112.md, from NOTE rw-p333 + rw-edd)
- About half of the reported omega_inf offset is a 1-mode fit bias. At P 3.33 tau 10, and only there among the 18 cells,
  the thermal-diffusion mode (-2.32i) decays slower than the acoustic mode (-3.66i). A 4-mode fit (acoustic, slow,
  left acoustic, fast; rms 1.6e-7) of viper's VET ladder gives the true acoustic omega_inf = 13.34589 - 3.66194i,
  off vetqs_d by (-4.8e-4, -1.31e-3).
- Exactly one continuous term reproduces the acoustic omega, the acoustic eigenvector (v, T, E, F phase -1.0e-4,
  -2.7e-4, -2.7e-4, -3.9e-4 rad vs rho) and both 1-mode limits (VET 13.34565 - 3.66301i, Eddington 13.70051 - 4.08921i,
  to 1e-5). That term is the momentum-exchange drag G = rho kappa_F F0/c, scaled by (1 + eps), eps = -3.15e-4, on
  BOTH sides (gas momentum and F equation). It is momentum-conserving and must shrink with P (consistent with
  eps proportional to E0/(rho c^2)). The analytic system is the exact linearisation, so the code is the side that is off.
- Suspect: implicit_gas_src = conservative2 + implicit_vimp (ImplicitCs2Vimp / Cs2DeltaA). The arms below localise it,
  and the c / P arms measure how eps scales.

## Build
FINAL 84cd20888 (the wt_rwa6_fin build of rw-amp1e6 / rw-edd; build_none, MPI, Release, -ffp-contract=off).

## Runner (docs/handover/scripts/rw_drag + rw_amp1e6 as siblings, as before)
rw_amp1e6/run_radwave.py now also writes the 4-mode acoustic omega to results.json (w4_re, w4_im, w4_slow,
w4_slow_amp, w4_eig_ratio; module rw_amp1e6/modefit4.py, needs scipy). It hashes tags over 250 characters (your
11-12 patch). Every case has ny=4 (nx2 = 4, as rw-edd).
```
cd <run area>/rw_drag
export OMP_NUM_THREADS=1 RW_PKF=$PWD/../rw_amp1e6/prodkeys_all.txt RW_EXE=<wt_rwa6_fin>/build_none/src/athena
RW_RUNS=$PWD/runs          python3 ../rw_amp1e6/drive.py run cases_drag.txt  --np 30 > cases_drag.out 2>&1
RW_RUNS=$PWD/runs_c20 RW_C=20 python3 ../rw_amp1e6/drive.py run cases_scale.txt --np 6 > cases_scale.out 2>&1
```
Pre-check: `python3 ../rw_amp1e6/drive.py input cases_drag.txt | grep -E "^time_scheme|^implicit_vimp =|^implicit_gas_src|dt_retry =|nx2"`
must show the arm's override, dt_retry = false and nx2 = 4. (rm_ keys come last in the input, so they override
the TIME production keys.)

| cases_drag.txt (N 64/128/256, nt = 2N, edd and vetf) | arm |
|---|---|
| vimp=true | base (control: must give eps ~ -3.15e-4) |
| + rm_implicit_gas_src=conservative | (i) the dA sharing of conservative2 off |
| rm_implicit_vimp=false rm_implicit_cs2_vimp=false rm_implicit_vimp_bc=false | (ii) explicit-velocity enthalpy flux |
| + rm_time_scheme=be | (iii) backward Euler (1st order in dt: eps from the N ladder is then dt-biased; compare its trend only) |
| P = 1.6667 | scaling in E0/(rho c^2): predicted eps x 1/2 |
| cases_scale.txt at RW_C=20 (P 3.33) | predicted eps x 1/4 |

A FATAL is a result (e.g. conservative may refuse a cs2 key): report it, do not add keys. Cost: N256 nt512 about 25 s per core.

## Return (NOTE-2026-11-12-orion-viper-rw-drag.md on this branch)
`python3 ana_drag.py runs cases_drag.txt` and `RW_C=20 python3 ana_drag.py runs_c20 cases_scale.txt`: paste both
whole (per N the 1-mode and 4-mode omega, slow-mode amplitude; the N->inf 4-mode limit, dw4 and eps(G)). Wall per list. NOTE only, no run trees.
