# TASK for Caltech: take over the accretor (Plaskett stream) project from DeltaAI

**User 2026-10-08:** the accretor project moves to Caltech. DeltaAI does no more accretor work. viper stopped on 10-07
(NOTE-2026-10-07-viper-accretor-handover). Nothing of this project is running or queued on DeltaAI.

## 1. Read first

1. `docs/handover/accretor-handover/INDEX.md`: viper's design, implementation and run history (RY Per to Plaskett).
2. `NOTE-2026-10-07-deltaai-plaskett.md`: q12_s1 / q12f_s1 (old setup) and their results.
3. **`NOTE-2026-10-07-deltaai-plaskett-env13.md`**: the setup audit, the bugs, the env13 fixes, validation, and the
   q13_s1 results. This is the current state.
4. `plaskett-1007/env13/RESULTS.md` + `plaskett_setup_1007.py`: the derivation of every env13 key.
5. Results page (private; ask the user for access if needed): https://claude.ai/artifact/RLmHQrvfFAr9yAv4eki4Tx

## 2. State

- **Physics setup:** env13 (`plaskett-1007/env13/plaskett_env13.athinput`, md5 313d7f0b1156519a8eaa7063162d07bc).
  - Envelope top on a Roche equipotential.
  - WB cut 0.41 Rsun below the photosphere on an equipotential.
  - Fluxes measured at r_meas 9.615675.
  - Ryu+2025 stream width 0.735 and env_rho_ph 0.055.
  - Wade+2026 temperatures.
  - 640 x 4 x 2048 grid.
  - **Spin 1.0, by user decision (10-07): no fast-spin arm.** 7.2 was RY Per's value, above Plaskett's critical 4.72.
- **Code:** branch **`accretor-1007`**, physics commit **`8cecb89a`**. Later commits are docs only. The new keys
  (env_top_mode, env_wb_depth, r_meas) are off by default and bitwise-neutral.
- **q13_s1 (DeltaAI, half orbit, done):**
  - Accreted fraction 0.919 (q12: 0.46).
  - No outflow through r_meas at any time.
  - Envelope edge steady to 0.013 Rsun.
  - j of accreted gas 0.556 j_K, which equals the ballistic orbit's (0.555).
  - 0 FATAL, min dt 9.58e-7.
- **Open (the user decides):** longer runs toward steady state. The input's own tlim is 10 orbits (4.58266).
  - Cost: q13_s1 ran at 2.66e8 zone-cycles/s on 4 GH200 (~51 cycles/s, 52 min per half orbit), so 10 orbits is about
    17 h on 4 GH200.
  - The pgen reads no files, so a fresh start needs only the input. The q13_s1 restarts (ryper.00000-00005.rst,
    about 400 MB each, the last at 0.5 orbit) and the bin dumps are on DeltaAI at
    /work/nvme/bivj/jma20/plaskett_1007/test1007/run/q13_s1/ if the user wants to continue rather than restart.

## 3. Do now

1. **Build 8cecb89a** for H200:
   `docs/handover/caltech-2026-09-26/scripts/build_inc.sh <tag> gpu 8cecb89a ry_per_accretor` (as in
   TASK-2026-10-07-plaskett-remote section 2). Record the md5. The DeltaAI GH200 md5 is dbbef84bb366998914264b111a8a22d3.
2. **Smoke:** env13, 20 cycles, 1 node x 4 GPUs (8 MeshBlocks, 2 per GPU), `-i plaskett_env13.athinput time/nlim=20`.
   The DeltaAI reference (job 3334844, t3_smoke) is in the table below. Expect agreement to ~1e-6 relative, not
   bitwise (FMA). A difference of 1e-3 or more is a port bug: stop and report.
3. **Report and wait.** Push `NOTE-2026-10-0x-caltech-accretor.md` to `accretor-1007` (fetch, merge, push; never
   force): md5, smoke numbers next to DeltaAI's, and cycles/s. **Do not start production until the user says what to
   run** (10-orbit q13 continuation or fresh start, layout, cadence).

DeltaAI reference values (job 3334844, t3_smoke):

| cycle | time | dt |
|---|---|---|
| 0 | 0 | 1.918640e-06 |
| 10 | 1.900746e-05 | 1.851290e-06 |
| 20 | 3.712850e-05 | 1.786984e-06 |

The run must report rc 0 and 0 FATAL. Startup lines to match:
- `measuring face (MR, JR, Menv, Jenv) at r = 9.61568`
- `envelope top on the equipotential psi = -15.0 c_ph^2 ... r_top 9.1635 at phi 90, max over phi 9.6863`
- `c_s,don 18.679 km/s ... width 0.7350 Rsun, sigma_phi 0.05444 rad`

## 4. Tools

- **Figures + page:** `plaskett-1007/env13/page/make_all.sh` with `mkfigs13.py` and `build_page13.py`.
  - Edit the paths at the top of the scripts.
  - Use the accretor-1007 `vis/python/bin_convert.py`: the rt-integration one cannot read bin format 1.2.
  - mkfigs needs numpy, scipy and matplotlib.
- **Numbers:** `plaskett-1007/env13/results/numbers_{q13_s1,q12_s1,q12f_s1}.json`.
