# Viper session state 2026-09-27 (checkpoint before a fresh session) — START HERE on viper

Read this, then `MEMORY.md` (its RULE ZERO points here), then `HANDOVER-2026-09-26.md` (the full
project state, shared with the Caltech session). Always `git fetch fork` before any push: the Caltech
session pushes to the same `rt-integration`. Never force-push.

## Running Slurm jobs (they survive the session)

| job | what | where | when done |
|---|---|---|---|
| 11995741 | WASP-121b **1x** continuation, rot 60 -> 300, 2 GPUs, 15 h | /viper/ptmp2/jinma/w121prod_0927/w1x | ~11 h after start |
| 11995742 | WASP-121b **10x** continuation, rot 60 -> 300, 2 GPUs, 23:55 | .../w10x | ~18 h after start |
| 11995660 | ke-dt: sph_wedge growth rate with the coupling fix | /viper/ptmp2/jinma/kedt_0926 | 1 h |

- Production binary: rt-integration 041fac8f, md5 7ea1537a (ck-jlin in, ck-store NOT in; ck-store is
  only -4 % at nx1 76).
- Run script: `run.sub` (restart from newest rst, STOP-file logic); the submit lines are in SUBMITTED.txt.
- At rot 60 the deep layer was still drifting (100 bar +13-16 K per 10 rot) -> target 300 rot on the
  low-res grid BEFORE any remap (user).
- 10x: 45 % of ck calls NOT-CONVERGED at rot 60 (max residual 8e-6, median 1e-7), with the 10x solver
  keys. Acceptable for relaxation; the denser-Jacobian fix is needed before a science stretch.
- Analysis: `w121prod_0927/ana_relax.py` (relaxation), `v60/` (verbose ck check).

## Open threads (each has a RESUME.md; the agents ended with the session)

| thread | branch | folder | state |
|---|---|---|---|
| **M1 coupling fix** (KE dt dependence) | ke-dt-0926 | /viper/ptmp2/jinma/kedt_0926 | gates (a) dt-independent growth 2.15e-4 at every cfl, (b) 2nd order in all variables: PASS with `force_reference_work=split` + `implicit_opac_update=true` + `implicit_one_pass=0` + `implicit_tol=1e-8`. Cost ~2.3x/cycle but ~15x cheaper at equal accuracy. Pending: sph_wedge (term + guard 81e88e44; job 11995660), LE shocks, tst, interleaved timing; branch tip 8f01828f, RESUME section "CHECKPOINT 09-27" -> then the merge/default decision (user) |
| **dhj flux history** | dhj-fluxhst (5f92ba50) | /viper/ptmp2/jinma/fluxhst_0927 (RESUME.md) | adds top/bottom radiative, energy and mass flux columns (`problem/flux_hst`); WIP, gates pending; merge when gated |
| **M1 + MHD coupling** | m1-mhd (52f1db4b, design only) | /viper/ptmp2/jinma/m1mhd_0927 | **ON HOLD (user 09-27) until ke-dt-0926 is merged**; then rebase on rt-integration and build on the merged fix. rad_m1 couples to hydro only today (no pmhd); Cartesian + sp; multi-day |
| remap nx1 76 -> 256 | — | radab_0924 (remap.py + remap_pgen.patch) | not started; needed only after rot 300 |
| 10x ck denser Jacobian | — | docs/handover/RESUME-ck-newton10x-b.md | not started |

## Done and merged this session (details in HANDOVER-2026-09-26.md)
- cs seam fixes, default on: positive / linear / guard 4 (root cause of all WASP-121b crashes).
- m1-wedge: `m1_test = sph_wedge`, top sponge default, zero-flux walls, open top outflow-only,
  `wg_bc_bot = reservoir`, super-Eddington study.
- ck-jlin and ck-store (from Caltech): viper HIP gates PASS.
- The Caltech port HIP gate PASS; kokkos submodule is now 4.6.02 (works for GFX942_APU).

## User rules added this session (full text in memory)
- Accuracy first, then speed; the accuracy bar is the ROUND-OFF spread, except the step size.
- Two phases: cheap schemes allowed in relaxation, accurate ones at steady state.
- Top layers don't count for accuracy (massive stars: tau < 1; dhj: p < 1e-6 bar), but they must not limit
  dt, crash, or cost much.
- Benchmark defaults on rad-hydro with moving gas, not static tests.
- Every fix: hydro AND MHD, Cartesian AND sp.
- Viper nodes have 2 GPUs; apu nodes are exclusive, so pair 1-GPU runs; size job limits ~1.5x the
  expected wall time.
