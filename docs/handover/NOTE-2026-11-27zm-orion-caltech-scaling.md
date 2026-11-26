# NOTE 2026-11-27zm orion -> Caltech: strong-scaling breakdown of a bsg cycle (rays vs VET solve vs hydro), USER ok

Question (USER): on many GPUs, does the ray sweep (2 x N_r serial shell steps, each with a lateral band exchange of ~90 us at 4 GPUs) come to dominate the VET solve (Krylov reductions + halos + mg coarse levels)?

Binary: release rel6 01deec542 (= aa2646b92 src), your existing build. Deck/state: bsg production restart (the rst you used in 27zc, cycle 29657), 20 cycles per arm, output off as in 27zc (time/ncycle_out=1 so elapsed lines are written), rad_vet/clean_verbose=1.
Arms: 1, 4, 8, 16 GPUs (1, 1, 2, 4 nodes), same total mesh, MeshBlock layouts that the release accepts (one MeshBlock along x1, whole tiles: same layout_tile in all arms so the answers stay comparable). If 16 does not divide the lateral mesh, take the nearest legal count and say which.
Report per arm (per cycle, mean over cycles 5-20):
- total s/cycle (elapsed lines) and zone-cycles/s;
- rays: the vet_gd 'fenced seconds' line split: source, sweeps, of which per-shell halo+walls, of which MPI wait; builds per cycle; shells swept (cut shell);
- VET solve: Krylov iterations per cycle, wall of the solve (clean stage steps), and if your kernel timer (kp) or nsys is at hand: time in Krylov kernels vs reductions (MPI_Allreduce) vs halo exchanges vs mg levels (per level if possible);
- hydro + rest = total - rays - solve.
One table + one plot (stacked bar per GPU count; mesa2.mplstyle + idort_style.py). Budget about 1 node-hour. Delete run dumps after.
