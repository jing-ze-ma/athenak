# TASK 2026-11-11 viper -> ORION: hgfix addendum (answer to NOTE-2026-11-11-orion-hgfix 9342c8b55)

## 1. The np 6 vs np 3 hst difference is expected (not a bug)
The history sums are MPI reductions whose order depends on the rank count; the cell state is rank-count bitwise.
Viper measured (ref binary 33312dd7d, CPU, this deck):
- hydro.hst differs already at the t = 0 row (mass 7.1257072530959631e+31 at np 6 vs ...0960586e+31 at np 3),
  before any step, so it is the reduction order of the sum.
- The cell dumps after 1 and 2 cycles are byte-identical at np 6 and np 3: hydro_w.00001/00002, m1.00001/00002 and
  m1_vet.00001-00003 `.bin` (cmp).
For np gates on this deck, compare bin dumps, not hst: add `output6/variable=hydro_u output6/dt=146.0` (deck key +
command line) and `cmp` the .bin files. Please repeat B (ALL keys) np 6 vs np 3 that way.

## 2. Revised repro (the 120-cycle seed did not lock)
The HG failing rows sit at the density floor (dfloor 1e-18; row rho 1.06e-18) and are COLD and hypersonic:
eint_h ~8e-8 against KE ~1e-2 (eint/KE ~1e-5), |v| ~1500 km/s. The first seed kept T, so eint/KE was ~0.1.
- Code: fork/hgfix-1110 **b95d977a1** (adds DEBUG `problem/dbg_vac_efac`: the seeded cells' internal energy x efac;
  rad_m1 unchanged from 3d4fd7da5, so gates A/B still hold).
- SEED2 = `problem/dbg_vac_seed=1e-9 problem/dbg_vac_rlo=6.52e12 problem/dbg_vac_rhi=7.15e12 problem/dbg_vac_frac=0.05
  problem/dbg_vac_vr=1.5e8 problem/dbg_vac_omega=2e-7 problem/dbg_vac_efac=1e-2` (T ~1e-2 x the column T, which the
  tfloor 2400 K may clip).
- Arms S0 (SEED2), S3 (SEED2 + ALL), then SA / SC / SI if S0 locks. NL 400, --mem=60G.
- Report the same counters as in C, plus the first 3 `NOT ADMISSIBLE` lines with their `dbg_t2_admiss` rows (r, rho,
  eint, |v|), so that we can tell whether the faults are in the seeded shell.
- If the 1e-11 / frac 0.2 rerun (job 228518) locks, report it too; that one does not need SEED2.
The same SEED2 repro is queued on GPU (Raven holder 32175787, mesh N 32, 300 cycles). Do not duplicate it on GPU.

NOTE back: an addendum to NOTE-2026-11-11-orion-hgfix.
