# NOTE 2026-11-25j orion -> Caltech: one stack capture of the nsys hang (after 25i)

Two of two nsys runs (mpi trace, then cuda only) stop in cycle 794785; HGN, HGO and HGV pass it. That may be nsys, or a wait in the code (the GPU ray halos with IPC, spin waits between ranks) that a slower timeline exposes, which we would want to know. After 25i is done, on 02-03 + 02-09:
1. nsys -t cuda again, same command, 3 cycles; when the log has not moved for 2 min: `gdb -p <pid> -batch -ex "thread apply all bt 12"` (or pstack / py-spy-free equivalent) of two ranks on each node (4 stacks), and `nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv` on both nodes; then cancel your own step by its id. If gdb attach is not permitted, say so and stop.
2. Same nodes, no nsys, same 3 cycles with `rad_vet/vet_gd_halo_ipc=false` is NOT asked; only (1).
3. If a Kokkos kernel-timer library (kokkos-tools kp_kernel_timer) is at hand for this CUDA build, 10 cycles of HGN's command with KOKKOS_TOOLS_LIBS set gives the kernel table without nsys; post it (text). If it is not there, do not build one.

-- orion (coordinator)
