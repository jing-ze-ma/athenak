# TASK raven -> Caltech: REAL BSG true reproduction with the FINAL method (bsg3d_truerepro2, 57.4 d), 1 node x 2 H200

From the Raven coordinator (user GO). Same run is queued on Raven (4 nodes x 4 A100, jobs 31107345-47, est. start ~20 h)
and on Viper (MI300A, 12164143-45). **First to start wins**: push a STARTED line in the NOTE as soon as link 1 runs; Raven
cancels the other copies (and tells you if another site started first: then cancel yours).

## 1. Build (report md5)
Fork branch **`bsgcand-1012` @ `d3bf42d7`** (= prodfinal-1011 ed32f0fa signed-off method + vet_gd_rotate_slope=correct +
jse face fix). Same toolchain/flags as your `_nofma` he_star_m1 H200 builds (`-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON
-D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D CMAKE_CXX_FLAGS=-ffp-contract=off`). Name `athena_gpu_he_star_m1_d3bf42d7_nofma`.
If the coordinator pushes a newer commit for this TASK before you build, use it.

## 2. Files
- IC + tables: your 10-09 BSG files dir (bsg-hrdet-1009 SETUP.sh; md5 ic_ma2026_f0 246a666e, rosseland 450bc0c1, planck e18ca805).
- Deck: `docs/handover/bsgprod-1012/bsgP_seed.athinput.in` -> `sed "s#@FILES@#<files dir>#g"`. It is the Raven production deck
  (original reproduction IC pinned, seed 1e-6 with he_seed_rad true, every final-method key named, rotate_dt 1.75e3 aligned,
  vet_gd_rotate_slope correct, implicit_thin_js_face hcface, mg_gf, tlim 4.96e6) with only: file paths templated,
  vet_gd_halo_ipc true and vet_gd_twin_lowmem true (single node, 8 blocks/GPU). Change no other key.
- `ic0check.py` + `ic0_ref_1009.json` (bsgrec-1011 package): cycle-0 fingerprint of the ORIGINAL IC.

## 3. Smoke, then production
- Smoke (`time/nlim=25 time/ndiag=1` - both in the deck): rc 0, no FATAL/NaN/NON-CONVERGED, dt 87.5 s every cycle, GPU peak
  (Viper used 85 GB/APU at 8 blocks), s/cycle; `python3 ic0check.py check <smoke dir> ic0_ref_1009.json`: dens/vel/m1_f* BITWISE,
  eint and m1_e may differ only in 37.6-45.1 Rsun at ~2e-6 / ~9e-6 (the seed, he_seed_rad = true; same as your 10-09 seeded run).
- Production: 1 node, 2 ranks x 8 MeshBlocks, 1 rank per H200 (your 10-09 layout), fresh t = 0 to tlim 4.96e6, chained links
  (restart from newest rst with rst_info.py last_time keys), STOP on rc != 0 / FATAL / NaN / NON-CONVERGED. Outputs as in the deck
  (hst 1e3 s, bins 0.5 d, rst 1 d). Expected ~1-1.5 s/cycle, ~57k cycles.

## 4. NOTE-2026-10-12-caltech-bsgprod.md (push to bsg-files-1009)
STARTED (job ids, binary md5, smoke numbers + ic0check output); after each link t, cycle, dt, s/cycle, Picard, counts, newest rst.
