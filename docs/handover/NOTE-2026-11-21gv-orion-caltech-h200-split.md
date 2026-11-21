# NOTE 2026-11-21gv orion -> caltech: 3 one-minute B_v203 runs to localise the H200 trap (while the rebuilds compile)

On a free test node (NOT 02-11/01-11 where the builds run), same tiny 1-GPU deck and B_v203 as the gq matrix, nlim 1.
Each run ~1 min. Report the clean_init pass 0-3 lines of each (trapped = pass 1 3.499e-08).
1. env M1_ARENA=0  (no memory arena; pass 1 is otherwise the first sweep into a rebound arena).
2. compute-sanitizer --tool initcheck  (uninitialised device reads; memcheck does not report these) — report the first
   20 errors with kernel names. If MPI_Init's UCX errors stop it, add UCX_TLS=^cuda_ipc,cuda_copy or run with
   OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader.
3. compute-sanitizer --tool racecheck  (shared-memory races), same handling; first 20 hazards with kernel names.
Nothing else; results to caltech_results/H200CHK/split.txt. The newer-CUDA rebuild (NOTE gt) stays the priority.

-- Orion coordinator
