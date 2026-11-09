# NOTE 2026-11-09 Orion -> viper: raycheap3 interim (the full NOTE-2026-11-09-orion-raycheap3.md follows)

Binary: 9fbd2d3c6 CPU build. src/, CMakeLists and kokkos are identical to 9b9951606. All arms x 4 levels were submitted at once (USER: "make tests fast").

**Done so far** (conv.py run by Orion on the finished levels; preliminary):
- **C1** (vet_col frozen: time2_vet_col = lag, vet_col_every = 1e6), all 4 levels:
  d 1.952e-4, 9.690e-5, 4.628e-5; **p 1.01, 1.07**.
  -> Freezing the closure keeps the scheme 1st order. The 1st-order term is NOT in the vet_col update. That points to the transport stepping or the forcing (your "next TASK" branch).
- **B0**, 3 coarse levels: d 1.805e-4, 8.932e-5; p 1.01 (reproduces viper). c0.01875 is at cycle ~2100 / 2400.
- **G1** (coupling = false), 3 coarse levels: d 2.529e-4, 1.878e-4; p 0.43. c0.01875 is at cycle ~2280 / 2400.

**C2 / C3 cannot run as specified.** Both stop at start, with or without `implicit_krylov_fuse=1`:
`rad_m1_sph.cpp:198: <rad_m1> on the cubed sphere ... supports only ... closure = eddington | vet_col ...; this input has: closure != eddington | vet_col`.
- So `closure = m1` is not available on the cubed sphere at 9b9951606. The fuse/jacobi workaround in the TASK does not apply.
- Substitute arm **C2e** (agent's choice): the C2 deck with `closure = eddington` + `implicit_hr_model = true`, as the nearest allowed no-VET baseline. Still running: c0.075 done; c0.0375 and c0.01875 in progress.
- Please say if you want a different substitute, or C2/C3 on the spherical-polar test star instead.

ETA for the full NOTE: ~30-60 min (C2e finest level).
