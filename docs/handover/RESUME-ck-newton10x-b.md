# RESUME: ck-newton10x-b (make the 10x ck Newton converge at tol 1e-8)

Branch `ck-newton10x-b`, worktree `/viper/ptmp2/jinma/wt_newton10xb`, based on rt-integration cab48eb7.
Background for all of this: `/viper/ptmp2/jinma/wasp121_0925/SMOKE_DIAG.md`, sections 4-7.

## State
- **Problem.** 10x WASP-121b with the production keys still has ck_implicit calls flagged NOT-CONVERGED.
  - Production keys: rsec 20, maxit 12, floorbound + kkt_demax, hiT2 tables, 1-D IC.
  - Flagged share: 16-19 % of calls at tol 1e-7 over 300 cycles (job 11980525), growing after cycle 250.
  - The slow columns are about 10 per call. They contract linearly, because the tridiagonal chord Jacobian drops the
    longer-range two-stream coupling. rsec only fixes the diagonal.
- **Commit 6702443b.** It removes the refusal of `ck_impl_aa` (Anderson/DIIS acceleration, already in the code,
  README_jac.md) together with ck_impl_floorbound / kkt_demax / rsec. That is the only code change so far.
- **GPU binary:** `/viper/ptmp2/jinma/bin_0925/athena.gpu.6702443b`, md5 4c8fb7e16d2d1e29216c6c2ac1c9d4c5.
  - Build script: `bin_0925/build_6702443b.sh`.
  - Never overwrite `sparc_0925/athena.gpu`: the 1x sponge runs use it.
- **CPU gates for 6702443b:** not run yet. They are needed only if the refusal change is kept.

## Option (b) design
1. **Anderson first (no new code):** `problem/ck_impl_aa = m` (depth 3/5/8), `ck_impl_aa_rst = true`, with rsec 20.
   - It needs `ck_impl_pred = false` (pred refuses aa) and `ck_impl_fuse = true`.
   - Anderson acts as GMRES preconditioned by the tridiagonal. It picks up the dropped coupling from the secant history
     at no extra sweep.
2. **If Anderson is not enough:** a denser Jacobian for flagged columns only, built from the stored factorisation
   (ck_impl_jac_lin / ck_lin_build keep lP/lG).
   - Either a banded or dense per-column two-stream response, or a per-column good-Broyden update on top of the chord
     tridiagonal.
   - Active only on columns still unconverged after pass 2-3.
   - New key, default off, bitwise.
3. **Keep defaults bitwise.** Accuracy may not change beyond the solver tolerance.

## Test job 11980859 (w121diag, apudev; submitted about 00:04)

| item | path or value |
|---|---|
| script | `/viper/ptmp2/jinma/wasp121_0925/smokediag/diag13.sub` |
| job log | `smokediag/log13.out.11980859` |
| arm dirs | `smokediag/<arm>/run.log` (+ bin, dhj.log, hst) |
| 10x input | `smokediag/in_b.athinput` = current `sparc_w121.athinput` (1-D IC, hiT2, rsec 20, maxit 12, floorbound + kkt) plus the aa keys |
| 1x input | `smokediag/in_b1x.athinput` = `x1/sparc_w121_1x.athinput` plus the aa / rsec keys |
| grid, 10x | `wasp121_0925/grid_w121.env` |
| grid, 1x | `wasp121_0925/x1/grid_w121_1x.env` |
| run settings | 400 fresh cycles, T4 + C2 + ck_impl_every 4, sponges on |
| common overrides | `problem/ck_impl_tol=1.0e-8` and `problem/ck_impl_debug=-2` (per-pass residual histories) |

Arms:

| arm | settings |
|---|---|
| b0 | production keys |
| b0np | b0 + pred off |
| aa3 | pred off + aa 3 |
| aa5 | pred off + aa 5 |
| aa8 | pred off + aa 8 |
| aa5m16 | aa 5 + maxit 16 |
| x1b0 | 1x, production keys |
| x1aa5 | 1x, pred off + aa 5 |

Analysis:

```
cd smokediag
python3 lev.py b0 b0np aa3 aa5 aa8 aa5m16 x1b0 x1aa5   # NOT-CONVERGED, residual median/max, passes, ms/cycle, histories
python3 ana.py <arms>                                   # floors, dt
```

**Pass gate:**
- 0 NOT-CONVERGED at tol 1e-8 over >= 400 cycles on the 10x arm.
- The 1x regression arm stays at 0 NOT-CONVERGED.
- Cost increase stays small (a few %).

## Exact next steps
1. Read job 11980859 with lev.py. Is any aa arm at 0/100 not converged at tol 1e-8, and is its convergence
   superlinear in the hist= histories? Compare ms/cycle with b0.
2. If an aa depth passes:
   - rerun it for >= 800 cycles, because the failures grew after cycle 250;
   - run the CPU gates (copy `/viper/ptmp2/jinma/newton10x_bld/gate/gate3.sh` pattern against a CPU build of cab48eb7:
     keys off bitwise, restart 6+6 bitwise with aa on);
   - commit, then report the recommended keys (aa m, pred off) for the 10x inputs.
3. If not: implement design step 2 (Broyden or a banded Jacobian on the flagged columns) behind a new default-off key.
   Then build into bin_0925, rerun the same arms plus the new key, and gate.
4. Update SMOKE_DIAG.md (new section 8) and commit. No merge, no push.

## 09-26: stalled columns diagnosed (job 11981384, binary bin_0925/athena.gpu.9326804a, md5 7fb751ee)
Script: `wasp121_0925/smokediag/diag14.sub`. Input: `in_s.athinput` (= in_b + the two new keys).
Runs: fresh 10x start, 200 cycles, 50 calls, tol 1e-8, maxit 12. Summaries with `python3 lev.py s_b0 ...`.

| arm | not converged / 50 | notes |
|---|---|---|
| s_b0 (production keys) | 14 | stalls, e.g. 4.89e-08 held for 4 passes with 8 columns left |
| s_nofb (floorbound off) | 13 | floorbound is NOT the cause |
| s_nokd (kkt_demax off) | 46 | residual pinned at 4.7e-5 |
| s_none (both off) | 46 | |
| s_krow (`ck_impl_kkt_row=true`) | 1 | the one left is cycle 0 (the start transient, 1.6e-7) |
| s_krowm16 (kkt_row + maxit 16) | 0 | 8.1 passes/call, 36.8 ms/cycle (b0 39.3; single runs, indicative) |

**Cause** (from the `ck_impl_stalldbg=6` dump, s_b0/run.log):
- Where: NIGHT-side columns, lon +-100..118 deg, lat +-21..70 deg (the dump shows mirror pairs).
- The stalled worst cell:
  - either i = 48-49 at ~1e-5 bar, T 230-250 K (just above tfloor);
  - or i = 52 at ~1e-4 bar, T 1280 K.
  - Its r/e is 1e-4 to 2e-3, i.e. about 3-5e-8 in the tol norm, and its applied step is 1e-11 to 1e-14 of e.
- The cell directly above it sits on the ck_impl_demax LOWER bound (kkt=1, cap=1). Its residual asks for cooling of
  0.4-1.8 e, and its applied step is 0.
- That KKT cell is removed from the convergence test, but its row still gives the tridiagonal an unclipped step of
  -0.4 to -1.7 e.
  - The neighbour's solution is built on that step through the off-diagonal.
  - The neighbour therefore converges to the fixed point of the clipped iteration, not to R = 0.
  - Its step decays quadratically to 0 while its residual stays constant: the stall.
- It is NOT round-off: hEm/e ~1e-4 and hS/e ~0.05 at the worst cell, so the round-off floor is ~1e-17 of e.
  Tol 1e-8 is achievable.

**Fix:** `problem/ck_impl_kkt_row` (commit f82d8591, namespace fix 9326804a; default false = bitwise).
- A thick KKT cell gets the identity row with d = 0, applied after rsec so the secant keeps the true residual.
  Its neighbours then solve the reduced (active-set) system.
- `problem/ck_impl_stalldbg = N` dumps every active column's worst cell, its neighbours and all clipped cells on passes
  >= N. Diagnosis only.
  - It is verbose: about 250 MB of log per 200 cycles, because a demax-clamped step differs from sx at round-off and
    so counts as clipped.

**Next:**
- Recommended 10x keys: `ck_impl_kkt_row = true`, maxit 16 (0/50 at tol 1e-8).
- Run the longer test (>= 800 cycles) plus the 1x regression arm.
- Run the CPU gates: keys off bitwise vs cab48eb7; 6 + 6 restart bitwise with kkt_row on.
- Both keys must be listed in the input before they can be overridden on the command line.
- The stalldbg run logs in smokediag/s_* are large (s_none ~1 GB) and can be deleted.

## 09-26 (2): gates and long arms (binary 9326804a; no merge)
**CPU gates** (`bin_0925/gate_9326804a`: gate.sh / gate.out, gate2.sh / gate2.out; CPU builds of cab48eb7 and
9326804a): all pass.
- Keys off vs cab48eb7: data arrays and hst are bitwise identical on sparc_w121 and on the old-planet sparc.athinput.
- Restart 12 cycles straight vs 6 + 6 with kkt_row on: bitwise.
  - In gate2 the cells are forced onto the bound with demax 0.02: kkt = 5600-6400 cells per call.
  - kkt_row on vs off differs there, which confirms the key is exercised.

**GPU, 10x, 800 fresh cycles, tol 1e-8, no stalldbg** (job 11981731; two repeats, interleaved side by side):

| arm | not converged / 200 | passes/call | ms/cycle |
|---|---|---|---|
| production | 164 | 11.27 | 43.3 |
| kkt_row + maxit 16 | 85 | 13.49 | 44.0 (+1.6 %) |

- Per 100 cycles, production fails 25/25 calls from cycle 200 on. kkt_row fails 0, 0, 5, 25, 19, 18, 10, 8.
- The failures left with kkt_row are NOT stalls. One or two columns are still contracting linearly at pass 16
  (1e-8 to 2e-7): the stiff day-side population. They need design step 2 (Broyden or a banded Jacobian).

**GPU, 1x:**
- Fresh start (job 11981732), 800 cycles: 9 -> 0 not converged, 34.25 vs 34.35 ms/cycle.
- From the production restart sparc_w121x1/base/rst/dhj.00017.rst (rotation 8.14, 1 rotation, job 11981784,
  `smokediag/rst1x`): 20/2113 -> 0/2112 not converged, 4.50 -> 4.34 passes/call, cpu time 287.9 vs 288.8 s.
