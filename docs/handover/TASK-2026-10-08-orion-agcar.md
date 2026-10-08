# TASK for Orion: try the AG Car 3-D shake-down on CPU (throughput test first, user 10-08)

## 1. Why
The user wants to know whether Orion CPUs are a usable third machine for the AG Car 3-D shake-downs (viper and Raven
GPU queues wait 5-15 h). Expected: one GPU node ~ 1000-3000 CPU cores for this physics, so measure before committing.
Copies of the same runs are queued on viper (12131308 B / 12131309 A) and Raven (31004544 B / 31004545 A).

## 2. Code and files
- Code: fork/rt-integration (67ab56ba or newer; contains he_ic_eint_from_t and everything AG Car needs).
  CPU build, PROBLEM=he_star_m1, MPI on, OpenMP on (Kokkos OpenMP backend), Release. Record md5.
- Files: this branch (he-ic-eint-from-t), `docs/handover/agcar-files-1008/`: ICs, TOPS tables, inputs with
  @AGCAR_DIR@, SETUP.sh, MD5SUMS, viper smoke refs in smoke_ref_viper/. Copy the directory to
  /orion/ptmp/jinma/agcar_1008/files and run `bash SETUP.sh /orion/ptmp/jinma/agcar_1008/files`.
- Use the inputs unchanged (they include rad_m1/implicit_opac_newton_slope_max = 3 and he_ic_eint_from_t = true).

## 3. Layout
The inputs have 4 MeshBlocks 480 x 64 x 64. nx1 MUST stay whole (one radial column per block: he_star_m1 / vet_gd).
For CPU you may split theta/phi: `meshblock/nx2=16 meshblock/nx3=16` (64 blocks) or 32x32 (16 blocks); the physics
is unchanged. Choose ranks x OpenMP threads per node yourself.

## 4. Smoke + throughput (B first, it is the short one)
1. `time/nlim=10` for B on 1 node, then on 4 nodes (same layout). Record zone-cycles/cpu_second from the log, s/cycle,
   and cores used.
2. Compare with viper smoke B (smoke_ref_viper/smokeB.hydro.hst, job 12130898): t, dt, mass (col 3), tot-E (col 7)
   after 10 cycles within ~1e-6 relative (CPU vs GPU is not bitwise), rc 0, no FATAL, no NON-CONVERGED, the two IC
   T-check lines ("IC column T(rho,eint)/T_col" ~6e-14, "he_ic_balance cells" ~1.08e-7).
3. Same for A on the larger of the two layouts (viper ref smokeA, job 12130897).
4. Project the wall time: B needs ~1,900 cycles at the starting dt (tlim 4.624e5 s, dt 242 s), A ~6,500 (tlim
   8.064e6 s, dt 1,253 s); the dt will fall once convection starts, so these are lower limits. viper GPU node: 2.5
   s/cycle; Raven node: 1.4-1.55 s/cycle.

## 5. Then STOP and report
Do NOT start the production runs yet. Write NOTE-2026-10-08-orion-agcar.md on this branch: binary md5, smoke
table vs viper, s/cycle and zone-cycles/s per core for each layout, node-hours projected for B and A, and how many
nodes you could realistically get. The user decides whether to run on Orion.
