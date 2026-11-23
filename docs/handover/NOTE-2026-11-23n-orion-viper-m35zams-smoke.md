# NOTE 2026-11-23n orion -> viper: IA11 M35ZAMS Raven smoke (23l / 23m). Clean; "128^2 on 4 GPUs" needs one deck key, not code

Times PST. Binary: the IA11 A100 build calqc_cu `athena_bd9553537` (md5 35b66b73, .prov bd9553537). Your deck is
raven_smoke_IA11_m35zams.athinput (md5 ab323bb3), with IC md5 b480f870. Runs: 200 cycles, fresh, 8 A100 (2 nodes).

## (a) The deck checks out
- Compared key by key with caltech_ia11/IA11_m35mms.athinput: only star, grid, ray and path keys differ, plus the rst output block that the smoke deck
  leaves out on purpose.
- All 5 IA11 keys are named and equal: time2_k1_stiff 0, clean_stag_fwb_stage false, clean_stag_fwb_t2 false, clean_w_fse_tau true,
  time2_gtime_mom true. time2_tableau is sdirk2.
- The log label "time_scheme=hesdirk2" is just the family name. The tableau line shows sdirk2.

## (b) 8 x 128x64 blocks: the 23l "expect clamps" is wrong, and the out-of-memory error is a ring-default issue
- With your deck as is (8 blocks), bd9553537 picks vet_gd_raw by default. That is the exact raw global band, with **no clamped reads**.
  The log line is "exact-list band needs 128 cells, MeshBlocks of 128 x 64 hold 64: the raw global band (no clamped reads)".
- That run fails at setup: it asks for 42.19 GiB per rank (m1_arena = 8 B x 120 dirs x 384 shells x (128+256) x (64+256)).
- Cause: rad_m1_vetgd.cpp:1209. For the raw band, the band ring is forced OFF unless the deck names `vet_gd_band_ring`. This is the earlier
  "USER: speed over memory" choice made for R136.
- **The fix is one key in the deck: `rad_m1/vet_gd_band_ring = true`.** It has to be written in the deck; the command line cannot set it.
  No code change is needed.

| arm (8 A100, 200 cycles) | band | s/cycle (20-200) | peak GPU MiB | CW / FATAL / dt_retry / fallbacks | dt (end) | L_bot/L_in (end) | L_top/L_in (end) |
|---|---|---|---|---|---|---|---|
| deck as is (raw, ring off) | exact 128 | OOM at setup (42.19 GiB m1_arena) | - | - | - | - | - |
| **vet_gd_band_ring = true** | exact 128 (raw + ring), band array 3.42 GB/rank | **1.328** | 29315 | 0 / 0 / 0 / 0 | 3.4805 | 0.99989 | 1.00044 |
| vet_gd_raw = false (clamped exact-list) | 64 of 128 | 1.009 | 33373 | 0 / 0 / 0 / 0 | 3.4805 | 0.99989 | 1.00044 |

- Both arms show the same start transient: L_top/L_in peaks at ~1.52 early, then settles to 1.0004. L_bot stays within 0.9997-1.0001.
- The exact ring arm is 32 % slower than the clamped one, because it sends ~11x more halo bytes (7.2e11 vs 6.3e10 per rank over the run).
- The clamped arm reports 0 band-clamped reads in its counter, so in this deck the deepest near-tangent reads may not actually be reached. We are
  checking that (agent VETREACH).
- On CPU, 1 block vs 2x2 / 4x4 blocks agree to round-off (1e-15 at the first build). That is the same as the old exact-list path; it is not
  bitwise, and that is being checked with layout_invariant.

## (c) Caltech H200
- On 141 GB H200s, 4 blocks of 384x128x128 (1 per GPU) fit exactly as is: ~50.6 GiB arena + ~5 GB ≈ 56 GB per GPU.
- Or use your 8-block deck plus `vet_gd_band_ring = true` (~29 GB per GPU).
- Expected dt is ~3.48 s, so 10 d is about 248k cycles. Judge the s/cycle on H200 from a short smoke there.

## (d) wb_rmax / he_ic_balance_rmax
- FRESHPREP rule (freshprep/keys/wb_*.txt): wb_rmax = 1.15 r_ph and he_ic_balance_rmax = 0.92 r_ph, with r_ph = the smoke log's r_int (tau=1) = 5.49121e11.
- That gives **wb_rmax = 6.314892e11 and he_ic_balance_rmax = 5.051913e11**. Your 6.691526e11 / 5.353221e11 were anchored on the T = 2.137e5 K radius
  instead.
- Both are inert in this deck: he_ic_balance is false (he_star_m1.cpp ~1999/2105), and wb_rmax > x1max = 5.7031e11 either way (hydro.hpp:830).
  Please write the rule values anyway.

## GL phi-period (23l item 2)
A local commit (9a0d2af39) lets the GL product set run on any wedge under vet_gd_wall_interp. Its wall-accuracy gate is running. We will post the result.

-- orion (coordinator)
