# TASK raven -> DeltaAI: REAL BSG true reproduction, final method, as a chain of <= 1 h links (1 node x 4 GH200)

From the Raven coordinator. **User GO for running this production as chained short links on DeltaAI** (exception to the
short-tests rule for this run). The same run is queued on Raven (31107345-47, est. ~20+ h), Viper (12164143-45) and Caltech
(TASK-2026-10-12-caltech-bsgprod.md). **First to start wins**: push a STARTED line as soon as link 1 runs; Raven then cancels
the other copies. If the coordinator posts that another site started first, stop after the current link.

## 1. Build (report md5)
Fork branch **`bsgcand-1012` @ `d3bf42d7`**, `build_inc_deltaai.sh hes_gpu` (PROBLEM=he_star_m1). `export MPICH_GPU_SUPPORT_ENABLED=1`.

## 2. Files
- IC + tables: your 10-09 BSG files (md5 ic_ma2026_f0 246a666e, rosseland 450bc0c1, planck e18ca805).
- Deck: `docs/handover/bsgprod-1012/bsgP_seed.athinput.in` (same as the Caltech TASK), `sed "s#@FILES@#<files dir>#g"`.
  It has vet_gd_halo_ipc true and vet_gd_twin_lowmem true (single node). Change no other key.
- `docs/handover/bsgprod-1012/ic0check.py` + `ic0_ref_1009.json` (cycle-0 fingerprint of the ORIGINAL IC).

## 3. Layout and links
- 1 node x 4 GH200, 4 ranks x 4 MeshBlocks (your bsgrec/bsgaccel N4 layout; ~50 GB/GPU with lowmem). Plain srun.
- Link 1: smoke-check the first 25 cycles inside the job (rc, no FATAL/NaN/NON-CONVERGED, dt 87.5 s, ic0check: dens/vel/m1_f*
  bitwise, eint/m1_e differ only in 37.6-45.1 Rsun at ~2e-6/~9e-6 from the seed with he_seed_rad = true), then keep running.
- Each link <= 1 h wall with `-t <wall - 5 min>` so the code writes a final rst; next link restarts from the newest rst with
  the rst_info.py last_time keys (bsg-hrdet-1009). Expected ~0.45 s/cycle -> ~57k cycles -> ~7-8 links to tlim 4.96e6 s.
  Submit links one after another as each ends (ghx4-interactive allows one job) or as a dependency chain on ghx4 if allowed.
- STOP rules: rc != 0 / FATAL / NaN / NON-CONVERGED -> stop, record cycle, t, last 40 lines of run.log; no rerun with changed keys.

## 4. NOTE-2026-10-12-deltaai-bsgprod.md (push to bsg-files-1009)
STARTED (job id, md5, smoke numbers, ic0check), then per link: t, cycle, dt, s/cycle, Picard mean/max, counts, newest rst.
Keep hst/rst/bins on DeltaAI; at the end post the run dir path and push user.hst + run logs.
