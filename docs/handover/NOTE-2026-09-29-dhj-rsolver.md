# NOTE 2026-09-29 (viper -> Caltech, DeltaAI): dhj uses lhllc / lhlld from now on; RESTART the Caltech WASP-121b runs

**User decision 09-29:** every deep_hot_jupiter_rt (dhj) run uses the low-Mach Riemann solvers from now on:
- `<hydro>/rsolver = lhllc` (was `hllc`)
- `<mhd>/rsolver = lhlld` (was `hlld`)

**Why:** the deep WASP-121b flows are low-Mach. On viper (mhdvh_0928) the plain HLL-type solvers damped low-Mach
flow 3-60x more than the low-Mach versions (He box: KE -67 % to -98 % with hlld vs lhlld), and the WASP-121b
production inputs were still on plain `hllc`.

## Caltech: restart the WASP-121b production NOW (user 09-29)

The runs started from NOTE-2026-09-28-caltech-w121prod.md (jobs 3607761 1x, 3607762 10x) use `rsolver = hllc`.
**User: restart both now, from scratch, with `rsolver = lhllc`.**
1. Stop both production jobs (and their chained links); keep the old run dirs for reference (rename, e.g. `_hllc`).
2. In both inputs (`/resnick/groups/carnegie_poc/jingze/w121prod_0928/in/w121prod_{1x,10x}.athinput`) set
   `rsolver = lhllc` in `<hydro>`. Nothing else changes (same commit/binary: `lhllc` needs no rebuild; Newton keys,
   xstep 8, closed wall, sponges, `flux_hst_floor` as before).
3. Smoke (50 cycles, fail on FATAL / NOT-CONVERGED) as before, then the fresh-start production chain to rot 300.
4. Push a one-line status to `docs/handover/` (or relay via the user): new job ids and the smoke result.

## Repo inputs switched (this commit)

- `docs/handover/caltech-2026-09-26/inputs/wasp121_{1x,10x}/*.athinput`: `lhllc`.
- `inputs/production/deep_hot_jupiter_cs_hyd4.athinput`: `lhllc`; `deep_hot_jupiter_cs_prod4.athinput`: `lhlld`.
- NOT changed: `docs/handover/bench-2026-09-28/inputs/*` (the cross-cluster timing benchmark is pinned; keep `hllc`
  there so all machines time the same thing).

## DeltaAI

Any dhj run you set up (beyond the pinned timing benchmark) uses `lhllc` / `lhlld`.

## Viper check (lhllc on WASP-121b, 400 cycles, fresh 1x and 10x, 2 GPUs, job 12018289)

RESULT: PASS. Both arms ran 400 cycles with rc 0, no FATAL, no NaN, 0 NOT-CONVERGED (with the Newton keys of the
bench inputs), Etot_bot/Lrad_bot = 1.000 and Mdot_bot = 0 (closed wall). The solver really changed: dt at cycle 20
differs from the hllc smoke in the 6th digit (1x 15.54181 vs 15.54189 s). Binary: 11c9a5be (md5 b8e3fff0); the
lhllc code path is identical at 6faba555 / e804daf1.
