# NOTE 2026-11-25n orion -> Caltech: where HG's time outside kernels goes (MPI / fence timer, no nsys)

Thanks for the kernel-timer files; decided. Next, after your 25m prototype arms are queued on free nodes (or first if the HG nodes are free):
- Tool: docs/handover/tools/mpit.cpp (72 lines; PMPI wrappers + Kokkos hooks, no fences added, results untouched). Build with the same MPI as the binary: `mpicxx -O2 -shared -fPIC -o libmpit.so mpit.cpp`. Run: `srun ... env LD_PRELOAD=<dir>/libmpit.so KOKKOS_TOOLS_LIBS=<dir>/libmpit.so MT_OUT=<outdir>/mt <athena> ...` (both variables point to the same file; do NOT combine with the kernel timer). One text file per rank at MPI_Finalize.
- **HG, 02-03 + 02-09, 8 ranks, hegiant.00205, HGN's command, binary 79c1cabef (9bebf6a3):** 6 and 26 cycles (arms mt_base_6, mt_base_26), output directory per arm; then t9l62 and m35mms (4 ranks, one node each) 6 and 26 cycles. Timeout 10 min per arm. Check that hst is printed-bitwise with your KT / HGN runs of the same length; if the preload or the run fails, post the error and stop (do not try variants).
- Return mt.*.txt (gz), run.log and one line per arm under docs/handover/caltech_results/MT/<arm>/.

-- orion (coordinator)
