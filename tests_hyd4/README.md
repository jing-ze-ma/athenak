# tests_hyd4 -- pure-hydro twin of the cs dhj prod4 input (2026-09-22)

Deliverable: `inputs/production/deep_hot_jupiter_cs_hyd4.athinput` (not committed, not launched).
Smoke/timing runs: `/viper/u2/jinma/ATHENAK/bench/cs_hyd4_smoke/` (submit.sh, runs/<arm>/run.log).
Binary: `/viper/u2/jinma/ATHENAK/bench/cs_mhd_prod4/athena` (md5 f95130b22e6b7f9f57d86ca44da22058,
395db5bc, HIP GFX942 + MPI, PROBLEM=deep_hot_jupiter_rt), used read-only.

## Does the setup run as pure hydro?

Yes, by construction: `MeshBlockPack::AddPhysics` builds Hydro when the input has `<hydro>` and
no `<mhd>`, and `src/pgen/deep_hot_jupiter_rt.cpp` branches on `phydro`/`pmhd` at every use
(boundaries, WB arrays, conduction, ck Rosseland table, sources, dumps); every `pmhd->` access
is behind a null check or the `is_mhd_` value.  WB (`wellbalance_dynamic`, `wb_x1`,
`wb_option`, `wb_cache_every`), `etotgrav`, all `rad_*` conduction keys and all EOS/floor keys
are read from whichever block exists (hydro.cpp:177-252, conduction.cpp, eos.cpp), with the
same defaults in hydro.cpp and mhd.cpp for every shared key (checked by diffing the two
GetOrAdd lists).

## Changes vs prod4 (H1-H6, also in the input header)

| | prod4 (MHD) | hyd4 | why |
|---|---|---|---|
| H1 | `<mhd>` | `<hydro>` | selects the module |
| H2 | rsolver hlld | hllc | hlld is MHD-only; hllc = the solver of sp_dhj_hyd and cs_hyd_rs/hllc; general EOS sets is_ideal (general_hyd.cpp:33) |
| H3 | cs_wellbalanced_src true | (absent) | **read only from `<mhd>`** (coordinates.cpp:42-49); under `<hydro>` it cannot be enabled with this binary, so hyd4 uses the cell-centre geometric source. It is exact for the isotropic pressure part and differs from the face-sum form only in rho v v. cs_hyd_rs/base.athinput had it under `<hydro>` where it was silently unread. The one place hyd4 is not "prod4 minus B"; enabling it needs a one-line change in coordinates.cpp. |
| H4 | ohmic_resistivity eos, max_eta 5e12, use_rkg_sts, cs_lowbeta_fallback | removed | MHD-only (resistivity.cpp and mhd.cpp read "mhd"; eos.cpp:296 would build the unused electron-fraction table if ohmic_resistivity were left in) |
| H5 | bbot 3.0, bc_outer_maxwell true | bbot 0.0, bc_outer_maxwell removed | bbot is REQUIRED for every hot_jupiter run (pgen.cpp:57/131 GetReal) and used only on MHD paths; removing it killed the first smoke job (11943300, "Parameter name 'bbot' not found"). bc_outer_maxwell only acts in the is_mhd_ branch. |
| H6 | output3 mhd_w_bcc | hydro_w | history becomes dhj.hydro.hst |

Kept unchanged: grid + stretch, spherical tm ck sweep + pseudo-spherical beam, WB polytropic /
wb_x1 / wb_cache_every 0, rot_potential (rotpot_src stays true), etotgrav, rad_implicit_x1,
rad_angular false, EOS table, pfloor 1e-5 barye / dfloor 5e-14 / tfloor_kelvin 200, dt_min
1e-3, cfl 0.3, rk2, all output cadences.  The hydro path has its own dt-collapse cell
diagnostic (Hydro::dt_diag, mesh.cpp:1043).

## Smoke test (job 11943303, apudev, 1 node, 2 ranks / 2 GPUs, 12 MB per rank)

Arms on one node (vipa1001), same binary, interleaved: M1 prod4 MHD 300 cycles -> HL hyd4 from
scratch, -t 00:07:40 -> M2 prod4 MHD 300 cycles -> H2 hyd4 300 cycles.  All rc 0, no FATAL.
(The first attempt, job 11943300, died at setup on the missing bbot; logs kept as
`log.out.11943300.bbotfail`, `runs_11943300_bbotfail/`.)

**It works.** HL ran 7000 cycles to t = 1.736e5 s (0.57 rotation) from scratch, wrote
rst 00000-00002 (incl. the 0.5-rot restart at 1.525e5) and hydro_w bins, and exited on the
wall-clock limit (runs/HL/run.log).

dt history (runs/*/run.log, ndiag lines):

| cycle | prod4 MHD (M1 = M2 = O12a) | hyd4 (HL) |
|---|---|---|
| 0 | 28.977 | 28.977 |
| 100 | 28.960 | 28.960 |
| 200 | 28.334 | 28.333 |
| 300 | 19.299 | 19.300 |
| 800 | -- | 17.61 (minimum of the run) |
| 1900 | -- | 23.61 |
| 3900 | -- | 25.20 |
| 6900 | -- | 26.03 |

The first 300 cycles match the MHD input to 4 digits (the 3 G field is not dt-limiting there);
afterwards hydro dt settles at 25-26 s (the MHD production is quoted at ~17-18 s, prod4 NOTES.md).

Event counters (runs/HL/dhj.log, per 3.05e4 s interval, ~1200 cycles each):

| cycle | dfloor | efloor | tfloor | vceil | eos_fail | fofc | tclamp |
|---|---|---|---|---|---|---|---|
| 1291 | 7.0e6 | 180238 | 1.84e6 | 0 | 0 | 0 | 0 |
| 2551 | 7.3e7 | 262304 | 6.15e6 | 0 | 0 | 0 | 0 |
| 3778 | 8.2e7 | 259964 | 6.04e6 | 0 | 0 | 0 | 8 |
| 5004 | 9.0e7 | 279695 | 5.92e6 | 0 | 0 | 0 | 0 |
| 6186 | 8.8e7 | 108823 | 5.46e6 | 0 | 0 | 0 | 0 |
| 7041 | 6.5e7 | 85278 | 3.91e6 | 0 | 0 | 0 | 12 |

No C2P failures.  dfloor is the same order as every earlier cs run's first intervals
(cs_hyd_rs/hllc 8.4e6 -> 9.4e7, cs_mhd_prod2 9.0e6 -> 7.6e7, cs_mhd_prod3 5e4 -> 4.4e7;
their dhj.log), i.e. the near-vacuum top.  tfloor (4-6e6 per interval) is ~20-60x those runs;
the likely cause is the pfloor 1e-3 -> 1e-5 change inherited from prod4, but there is no MHD
prod4 dhj.log beyond 300 cycles yet to compare with (the production chain 11941995-8 is
still pending), so this is NOT isolated.  efloor 1e5-3e5 per interval (prod4's CPU
reproducer predicted efloor ~0 for MHD with pfloor 1e-5; unverified for hydro).
Warnings: the Rgas-unused notice (also in the MHD arms), and in HL one
"explicit radiative source was clipped in 2 cell(s) by rt_de_max" between cycles 300 and 400.
History (dhj.hydro.hst): mass 3.60051e26 -> 3.59550e26 (-0.14%), total E 6.6463e38 ->
6.6272e38 (-0.29%) by t = 1.747e5.

Timing (GPU, same node, same binary; cycles/s from the elapsed= lines):

| arm | cycles 100-300 | cycles 1900-6900 |
|---|---|---|
| M1 prod4 MHD | 14.23 | -- |
| HL hyd4 | 15.97 | 15.41 |
| M2 prod4 MHD | 14.21 | -- |
| H2 hyd4 | 15.56 | -- |

Hydro is ~10% faster per cycle (15.6-16.0 vs 14.2) and takes ~1.4x larger steps after the
spin-up (25-26 s vs the ~17-18 s quoted for MHD), so ~1.5x more simulated time per GPU-hour.

## What a hydro production would need (not launched)

* A bench dir like bench/cs_mhd_prod4 (same binary, or a fresh HEAD snapshot; nothing
  in src/ that this run uses has changed since 395db5bc), with its submit.sh/chain.sh copied
  and the job name changed; input = this file.
* Cost at 15.4 cycles/s and dt ~25.5 s: 8.64e7 s = ~3.4e6 cycles = ~61 h, i.e. 3 chained
  24 h apu links (4 to be safe if dt drops as the atmosphere spins up).
* Output cadence as prod4: hst/log every 0.1 rot, hydro_w bin every 2 rot (16 MB each),
  rst every 0.5 rot (78 MB each; 566 over 283 rot = ~44 GB unless old restarts are pruned
  -- the same issue as prod4).
* Decide H3 first: to be a clean "prod4 minus B" control the cubed-sphere face-sum source
  (cs_wellbalanced_src) should be readable from `<hydro>` (coordinates.cpp:42-58); as it
  stands hyd4 uses the cell-centre geometric source.
* Watch tfloor against prod4's own dhj.log once that chain has run.
