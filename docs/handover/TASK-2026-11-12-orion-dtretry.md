# TASK 2026-11-12 viper -> Orion: dt_retry fix, CPU correctness gate (HG small cs deck)

## Why
The amvecH gate found a forced dt retry with View poisoning crashing (segfault in RadiationM1::T2AdmissDebug).
Viper diagnosed it on GPU (apudev, HG small deck 704x32x32, base 84cd20888) and reproduced it once on CPU:
- The poisoned View is m1.iw, and only cells that no kernel ever writes: the x1 ghost layer of M1_IW_EP (slot 1)
  and of the implicit_vimp x1 coefficients P1m/P10/P1p (iw_vimp+0..2). Their values are 0 before the poison (0
  nonzero cells, max |value| 0), so they hold 0 in every attempt and in every restart. The whole-View poison
  NaN'd them, and an all-NaN E made T2AdmissDebug's MinLoc return the identity location (out-of-range read). So
  this is a false positive of the test, not a stale value.
- The real stale state: under the TIME keys (implicit_hr_cmem + physflux), a forced retry at fac 1 was NOT
  bitwise equal to the straight run (base, GPU), because fax1..3 (the applied face fluxes, state across steps) were
  never visited by the dt_retry checkpoint. The generated visitors were also stale (t2_abe_*, impl_cm_fainit_done,
  xe_init_pending, rkg_lmax_done, dt_cell, ...).
- Fix branch **fork/dtretry-fix-1111 3ec401877** (off final-1111 84cd20888): regenerated visitors; fax1..3 added
  to the minimal list; dtr_ncut set after the put-back; T2AdmissDebug made NaN-safe; dt_retry_poison = 1 now
  poisons only the elements the failed attempt changed (= 2 is the old whole-View poison); debug keys
  dt_retry_poison_n / dt_retry_poison_cells.
- GPU (viper) results with the fix: HG deck forced retry == straight, bitwise, with the production keys and with the
  TIME keys; poison=1 runs clean and is bitwise equal to straight; keys absent bitwise equal to base. This TASK is the
  CPU counterpart (the viper CPU queue is jammed).

## Build (two binaries, same recipe)
```
git worktree add wt_dtr_fix fork/dtretry-fix-1111     # 3ec401877
git worktree add wt_dtr_base 84cd20888                # fork/final-1111
# in each:
cmake -B build_hf -D Athena_ENABLE_MPI=ON -D PROBLEM=he_star_m1 -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build_hf -j32
```

## Run (docs/handover/scripts/dtretry_1111; tables = docs/handover/hgpos-1109/files, MD5SUMS there)
```
cd <run area>; cp -r <repo>/docs/handover/scripts/dtretry_1111 . && cd dtretry_1111
BASE=<wt_dtr_base>/build_hf/src/athena FIX=<wt_dtr_fix>/build_hf/src/athena \
TB=<repo>/docs/handover/hgpos-1109/files ./run.sh > gate.out 2>&1
```
This runs 7 arms at once, 3 MPI ranks each (21 cores), 4 cycles, with a restart dump at cycle 4. Expected wall time
is about 30-45 min (on viper one 3-rank arm takes ~7 min per cycle).

| arm | binary | keys | expected |
|---|---|---|---|
| BN | base | production keys (prodkeys.txt) + seed | reference |
| CN | fix | same | == BN bitwise (keys absent) |
| CF | fix | + forced retry at cycle 3, fac 1.0 (views minimal = production) | == CN bitwise |
| CP | fix | + views=all, audit, poison=1 | 0 FATAL, no POS_TRIGGER, == CN bitwise |
| CTN | fix | + TIME keys (cmem, cmem_cap 1, x2_fields 15, physflux mom fa zold fc, cmem_rst, fainit) | reference |
| CTF | fix | CTN + forced retry | == CTN bitwise |
| CTP | fix | CTF + views=all, audit, poison=1 | 0 FATAL, == CTN bitwise |

## Return (NOTE-2026-11-1x-orion-viper-dtretry.md on this branch)
Paste the whole gate.out: the 7 arm lines (rc, FATAL, POS_TRIGGER, the poison line) and the 5 cmp blocks
(IDENTICAL/DIFFER per restart file). For any DIFFER, also paste `grep -E "DT_RETRY|dt_retry (poison|audit)" runs/<arm>/run.log | cut -c1-400`.
Commit the NOTE only (no run trees, no restart files: each is a few GB).
