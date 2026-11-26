# NOTE 2026-11-27zt orion -> Caltech: VET-solve sync diet, bsg s/cycle and GPU idle before/after (VETPAR)

From agent VETPAR (orion). Priority per 27zo: after 27zm/27zn. This is a timing measurement: the arms keep a node to themselves. Short runs only.

## Binaries
- **REL**: your athena_rel6_01deec542 (= aa2646b92 src).
- **VP**: branch **`caltech-vetpar-aa6b0e57d`** (pushed with this NOTE) = aa2646b92 + 6 commits (`git log --oneline aa2646b92..aa6b0e57d`):
  - A1: max-norm test inside the sums' reduction;
  - A2: merged count reductions;
  - A4: stream-ordered zero fills;
  - A5: DevCopy cache;
  - B1a/b: halo exchanges at depth 1.

  All six are meant to be bitwise.
- Orion CPU gates (jobs 243315, 243322):
  - builds OK (star_vet md5 354f2ea0, built_in_pgens c6fe7b22); tiny_bsg_rel5 3.720e-11 PASS for every commit;
  - 49 output files (dumps, hst, rst) byte-identical to aa2646b92 on tiny_bsg_rel5 (256x16x16, 4 MeshBlocks), star_small, shell_sp and shell_cs at 2 and 4 ranks, for cbd6d417f (A1-A5), 46137bb76 (B1a) and aa6b0e57d (B1b);
  - restart 4 -> 8 bitwise on shell_sp and shell_cs at 4 ranks;
  - quick suite 25/25.
- Build as REL: CUDA 12.9 + g++13, star_vet; never CUDA 12.2. Run the tiny gate first (3.720e-11) and post md5 + prov.

## Arms (bsg production restart, the one used in 27zc/27zm; outputs off as in 27zm; clean_verbose 1)
| arm | binary | GPUs | cycles |
|---|---|---|---|
| R4 | REL | 4 (1 node) | 30 |
| V4 | VP | 4 (same node type, same layout) | 30 |
| R4n / V4n | REL / VP under `nsys profile` (cuda, mpi, osrt), cycles 11-20 | 4 | 20 |
| R8 / V8, R16 / V16 | only if 27zm ran 8/16 GPUs: the same MeshBlock layout and `layout_tile` as its arms | 8 / 16 | 20 |

## Report
1. s/cycle (mean over cycles 6-30) for R4 and V4. Also the 8/16 pairs if run.
2. **Bitwise:**
   - the CLEAN lines (`it K (sum T) ... solves S`) of R4 and V4 must be identical, cycle by cycle;
   - the final hydro_w/vet dump of a 10-cycle pair (write one at the end) must be byte-identical. Say yes or no.
3. From nsys (rank 0, cycles 11-20), per cycle, R vs V:
   - GPU idle ms and the number of idle gaps;
   - counts of `cudaStreamSynchronize`, `cudaDeviceSynchronize`, `cudaMemcpyAsync` H2D and `MPI_Allreduce` / `MPI_Allgather` / `MPI_Waitall`;
   - Isend bytes.

   Expected for V vs R (4 GPUs):
   - Allreduce down by about 60-100 per cycle;
   - device syncs down by about 60;
   - H2D copies from about 1250 to a few hundred;
   - rad-vet Isend bytes down by about 30-50 %;
   - idle down by about 10-20 ms/cycle.
4. Delete run dumps afterwards. Budget: about 1 node-hour at 4 GPUs, plus the 8/16 arms if run.

## Addendum (VETPAR, same day): preconditioner-setup share (USER question)
Add **one arm R4k**: REL, 4 GPUs, cycles 1-20 of the same restart, under SCSPEED's fenced Kokkos kernel timer (`KOKKOS_TOOLS_LIBS=libkp_ktime.so`, as in 27ze/27zi).

Report per cycle, mean of cycles 6-20, rank 0, ms and launch count for these labels:
- `cln_fwc` (face weights, 3 per setup) and `cln_stencil` (the operator's 19 coefficients). Both are the operator itself.
- `cln_pcsetup` (the x1-line and lateral preconditioner rows) and `cls_galerkin` (the coarse levels). These are the only reusable part.
- For scale: `cln_pcr2`, `cls_restrict`, `cls_prolong`, `cln_apply_st` and `cln_apply_st_sh`, `vet_lay_red2`.
- The CLEAN lines' `it K (sum T) ... solves S` per stage.

These decide whether reusing the preconditioner across passes, across the two stages, or until K rises can pay off. A reused setup saves `cln_pcsetup` + `cls_galerkin`, and costs about 3 ms of GPU time per extra Krylov iteration.
