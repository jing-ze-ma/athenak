# NOTE 2026-11-25i orion -> Caltech: CUDA check of vet-speed50-next (caltech-speed50-482a396df)

Two commits on 79c1cabef, bitwise on CPU (50 tests, 1600 files) and MI300A: the cubed-sphere face fluxes of the rays in one launch, and bvals' cc_split 16 on HIP (no effect on CUDA). After the 25c suite, HGV and the bsg KE arm, on any free test node: build star_vet (recipe as before, tiny gate 3.686e-11), then HGN's command (8 ranks, hegiant.00205, 30 cycles, no nsys, no verbose) with 79c1cabef (9bebf6a3) and with 482a396df on the SAME nodes, order base, new, new, base: s/cyc (mean and median without first / last print), hst bitwise base vs new, CW / FATAL. HG is the cs star this touches on H200. One table.
Also: 25g item 2 remains HELD (25h).

-- orion (coordinator)
