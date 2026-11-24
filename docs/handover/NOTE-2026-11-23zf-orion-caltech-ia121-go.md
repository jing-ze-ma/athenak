# NOTE 2026-11-23zf orion -> Caltech: GO IA12.1 (USER-approved): rebuild on fa59bffd6, restart t9l62 from 00012, then the M35ZAMS smoke

The USER approved IA12.1 for Caltech. The final tip is `caltech-ia121-fa59bffd6` on GitHub = your 21e176c67 + one commit (5 lines in rad_m1_vetgd.cpp: lets GL 7x18 run on a narrow wedge; only a startup check changes).

## 1. Rebuild
- Incremental rebuild in place, as in 23ze (CUDA 12.9 + g++ 13).
- Tiny BSG gate: expect pass 1 = 3.686e-11 again. Stop and post if it differs.

## 2. t9l62 restart (hpc-sm-02-09, 4 H200)
- `athena -r t9l62.00012.rst -d <new dir> $(grep -v '^#' IA121_t9l62_restart.args | tr '\n' ' ')`, plus your usual output / time keys.
  - The args file is `docs/handover/caltech_ia121/IA121_t9l62_restart.args` (md5 in MD5SUMS). It names every IA12 / IA12.1 key, because a restart otherwise takes the rst's stored values.
- Use 00012 (t ≈ 4.752e5). NOT 00013 (after the NaN).
- New output dir. Keep rst_keep's hard links.
- Watch and post:
  - the first 300 cycles: CLEAN WARNING / FALLBACK / fallbacks / dt_retry, dt, s/cycle, Ltop/Lin, Lbot/Lin;
  - the `CLEAN WBFDIAG` and `floor_audit` lines through the outburst phase (t 4.75-4.9e5, where IA11 failed at 4.85e5).
  - Then the 2-h status line as before.
  - Alarm, stop and post: dt < 1, dt_retry > 3/h, any NaN / inf, or unpaid floor energy.
- **Scaling later:** t9l62 has 8 MeshBlocks with tile = block, so it can move from 4 to 8 H200 at any restart with an identical answer. When bsg finishes and its node frees, tell us. The USER decides whether t9l62 takes it.

## 3. M35ZAMS 128² smoke (23zc step 3)
- On the fixed build, once t9l62 is running, wherever 4 H200 are free (bsg's node when it finishes), as prepared in `smoke_ia121_1123/run_smoke.sh`.
- Note: each rst is ~9 GB (GL intensities), so KEEPRST=3 is ~27 GB.
- Our Raven A100 run of the same deck (2x2 of 64²) gives 0.57 s/cycle, 0 CW, L_top/L_in 1.00046 at cycle 199, and dt 3.48 then falling to 2.48 by cycle 300.

bsg and m35mms: keep running on IA11.

-- orion (coordinator)
