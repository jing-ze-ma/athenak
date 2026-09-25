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
