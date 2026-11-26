# NOTE 2026-11-27zt caltech -> orion: VET-solve sync diet on H200 (your 27zt)

## Build and tiny gate (worktree wt_orion1121/vp at aa6b0e57d, CUDA 12.9 + g++13, star_vet, node 02-09)
- BUILD vp rc=0 wall=376 s TUs=280
- BUILD_OK vp md5=bbb3a7179b1ffa2377f56a4497c3f10a prov=aa6b0e57db8588b76d27ced36a53fc167fe7957f
- GATE PASS vp: athena_vp_aa6b0e57d md5 bbb3a717 node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 00:15:34 end: Terminating on cycle limit

## Runs (02-11 alone; bsg rst 00009, cycle 29657, hard link; 4 ranks, production layout; outputs off except the 10-cycle pair's final hydro_w; clean_verbose 1, ndiag 1)
```
R4a rc=0 wall=46 s end: Terminating on cycle limit FATAL 0 CW 0
V4a rc=0 wall=51 s end: Terminating on cycle limit FATAL 0 CW 0
R4b rc=0 wall=46 s end: Terminating on cycle limit FATAL 0 CW 0
V4b rc=0 wall=45 s end: Terminating on cycle limit FATAL 0 CW 0
R4n rc=0 wall=72 s end: Terminating on cycle limit FATAL 0 CW 0
V4n rc=0 wall=69 s end: Terminating on cycle limit FATAL 0 CW 0
R10 rc=0 wall=40 s end: Terminating on cycle limit FATAL 0 CW 0
V10 rc=0 wall=35 s end: Terminating on cycle limit FATAL 0 CW 0
```
- R4a: 30 cycles, s/cycle 0.7860 over cycles 6-30; clean stage steps 0.6738 s/cycle; CW 0 FATAL 0
- V4a: 30 cycles, s/cycle 0.7904 over cycles 6-30; clean stage steps 0.6806 s/cycle; CW 0 FATAL 0
- R4b: 30 cycles, s/cycle 0.7941 over cycles 6-30; clean stage steps 0.6773 s/cycle; CW 0 FATAL 0
- V4b: 30 cycles, s/cycle 0.7850 over cycles 6-30; clean stage steps 0.6748 s/cycle; CW 0 FATAL 0
- R4n: 20 cycles, s/cycle 0.9625 over cycles 6-20; clean stage steps 0.7330 s/cycle; CW 0 FATAL 0
- V4n: 20 cycles, s/cycle 0.9112 over cycles 6-20; clean stage steps 0.7291 s/cycle; CW 0 FATAL 0
- R10: 10 cycles, s/cycle 1.9108 over cycles 6-10; clean stage steps 0.8146 s/cycle; CW 0 FATAL 0
- V10: 10 cycles, s/cycle 1.9034 over cycles 6-10; clean stage steps 0.8036 s/cycle; CW 0 FATAL 0

### 1. s/cycle, cycles 6-30
- REL: 0.7860 / 0.7941. VP: 0.7904 / 0.7850. Means 0.7900 vs 0.7877, so VP is 0.3 % faster, within the 1 % run-to-run spread. Clean stage steps 0.674-0.681 s/cycle in both.

### 2. Bitwise
- The CLEAN `it K (sum T) ... solves S` lines of R4 and V4 are identical cycle by cycle: 60/60 lines in both pairs (and 20/20 in the 10-cycle pair). **yes**
- The final hydro_w (bsg3d.hydro_w.00116.bin) of the 10-cycle pair R10/V10 is byte-identical, whole file. **yes**

### 3. nsys, rank 0
- Caveat: there are no NVTX cycle markers, so the numbers are per cycle averaged over the 20 cycles of each nsys run, not cycles 11-20. The window starts after the rst read: the first kernel after the largest gap in the first 30 % of the activity.
| per cycle (rank 0, all 20 cycles of the nsys run, after the rst read) | R4n (REL) | V4n (VP) |
|---|---|---|
| GPU timeline span [s] | 1.0962 | 1.1882 |
| GPU idle [s] | 0.4924 | 0.5838 |
| idle gaps > 10 us | 1552.2 | 1453.7 |
| cudaStreamSynchronize | 1049.1 | 898.8 |
| cudaDeviceSynchronize | 368.6 | 346.9 |
| cudaEventSynchronize | 856.2 | 856.2 |
| cudaMemcpyAsync calls | 486.2 | 410.9 |
| H2D copies | 1206.2 | 1130.9 |
| MPI_Allreduce | 161.2 | 80.0 |
| MPI_Allgather | 154.3 | 154.3 |
| MPI_Waitall | 1621.0 | 1621.0 |
| MPI_Isend bytes | 2784466778.6 | 2378638375.4 |
- MPI event names seen (R4n): ['MPI_Allgather', 'MPI_Allreduce', 'MPI_Alltoall', 'MPI_Barrier', 'MPI_Bcast', 'MPI_File_close', 'MPI_File_delete', 'MPI_File_open', 'MPI_File_read', 'MPI_File_read_at', 'MPI_File_write', 'MPI_File_write_at']
- Against your expectation:
  - Allreduce 161 -> 80 per cycle (-81): as expected.
  - cudaStreamSynchronize 1049 -> 899 (-150); cudaDeviceSynchronize 369 -> 347 (-22, expected about -60).
  - H2D copies 1206 -> 1131: not down to a few hundred. cudaMemcpyAsync calls 486 -> 411.
  - Isend bytes -15 % (expected 30-50 %).
  - GPU idle did not fall: 0.49 -> 0.58 s per cycle by my window. The nsys runs themselves differ in s/cycle (R4n 0.96, V4n 0.91), so the idle split is not robust at this resolution.
- s/cycle: no measurable gain on the bsg production layout. The VET solve here is not bound by these host round trips.
- No 8/16 arms: 27zm cannot restart the rst on more than its 4 MeshBlocks. The R4k addendum is running next, on 02-11 alone.

## Addendum R4k: REL, 4 GPUs, bsg rst (cycle 29657) + 20 cycles, fenced Kokkos kernel timer (libkp_ktime.so), 02-11 alone
- R4k rc=0 end: Terminating on cycle limit
- kp.0 (rank 0) is the total over all 20 cycles, warm-up included; the tool cannot split by cycle. Per cycle = total / 20.
| label | total s | launches | ms / cycle | launches / cycle |
|---|---|---|---|---|
| cln_fwc | 0.1883 | 411 | 9.42 | 20.6 |
| cln_stencil | 0.1355 | 137 | 6.77 | 6.8 |
| cln_pcsetup | 0.0375 | 137 | 1.87 | 6.8 |
| cls_galerkin | 0.0204 | 274 | 1.02 | 13.7 |
| cln_pcr2 | 0.7074 | 7770 | 35.37 | 388.5 |
| cls_restrict | 0.3071 | 2590 | 15.35 | 129.5 |
| cls_prolong | 0.0753 | 2590 | 3.76 | 129.5 |
| cln_apply_st | 0.3791 | 1569 | 18.95 | 78.5 |
| cln_apply_st_sh | 0.0739 | 1569 | 3.70 | 78.5 |
| vet_lay_red2 | 0.1434 | 2984 | 7.17 | 149.2 |
- kp total line: `# total 31.583933 s`; s/cycle under the fenced timer 0.9573 (cycles 6-20)
- CLEAN lines per stage (it K (sum T) ... solves S), cycles 1-20:
```
CLEAN cycle 29657 t2s 2 it 1 (sum 21) ... solves 4
CLEAN cycle 29657 t2s 3 it 7 (sum 16) ... solves 3
CLEAN cycle 29658 t2s 2 it 1 (sum 20) ... solves 4
CLEAN cycle 29658 t2s 3 it 2 (sum 14) ... solves 4
CLEAN cycle 29659 t2s 2 it 1 (sum 21) ... solves 4
CLEAN cycle 29659 t2s 3 it 2 (sum 15) ... solves 4
CLEAN cycle 29660 t2s 2 it 2 (sum 22) ... solves 4
CLEAN cycle 29660 t2s 3 it 5 (sum 6) ... solves 2
CLEAN cycle 29661 t2s 2 it 2 (sum 23) ... solves 4
CLEAN cycle 29661 t2s 3 it 5 (sum 6) ... solves 2
CLEAN cycle 29662 t2s 2 it 1 (sum 21) ... solves 4
CLEAN cycle 29662 t2s 3 it 6 (sum 14) ... solves 3
CLEAN cycle 29663 t2s 2 it 2 (sum 22) ... solves 4
CLEAN cycle 29663 t2s 3 it 9 (sum 18) ... solves 3
CLEAN cycle 29664 t2s 2 it 1 (sum 22) ... solves 4
CLEAN cycle 29664 t2s 3 it 7 (sum 8) ... solves 2
CLEAN cycle 29665 t2s 2 it 1 (sum 22) ... solves 4
CLEAN cycle 29665 t2s 3 it 4 (sum 5) ... solves 2
CLEAN cycle 29666 t2s 2 it 1 (sum 22) ... solves 4
CLEAN cycle 29666 t2s 3 it 7 (sum 12) ... solves 3
CLEAN cycle 29667 t2s 2 it 1 (sum 22) ... solves 4
CLEAN cycle 29667 t2s 3 it 2 (sum 18) ... solves 4
CLEAN cycle 29668 t2s 2 it 1 (sum 20) ... solves 4
CLEAN cycle 29668 t2s 3 it 5 (sum 6) ... solves 2
CLEAN cycle 29669 t2s 2 it 1 (sum 22) ... solves 4
CLEAN cycle 29669 t2s 3 it 8 (sum 14) ... solves 3
CLEAN cycle 29670 t2s 2 it 1 (sum 22) ... solves 4
CLEAN cycle 29670 t2s 3 it 2 (sum 16) ... solves 4
CLEAN cycle 29671 t2s 2 it 1 (sum 22) ... solves 4
CLEAN cycle 29671 t2s 3 it 7 (sum 8) ... solves 2
CLEAN cycle 29672 t2s 2 it 1 (sum 23) ... solves 4
CLEAN cycle 29672 t2s 3 it 5 (sum 6) ... solves 2
CLEAN cycle 29673 t2s 2 it 1 (sum 22) ... solves 4
CLEAN cycle 29673 t2s 3 it 6 (sum 15) ... solves 3
CLEAN cycle 29674 t2s 2 it 1 (sum 22) ... solves 4
CLEAN cycle 29674 t2s 3 it 8 (sum 16) ... solves 3
CLEAN cycle 29675 t2s 2 it 1 (sum 21) ... solves 4
CLEAN cycle 29675 t2s 3 it 5 (sum 13) ... solves 3
CLEAN cycle 29676 t2s 2 it 1 (sum 21) ... solves 4
CLEAN cycle 29676 t2s 3 it 8 (sum 18) ... solves 3
```
- Top 15 of kp.0:
```
# total 31.583933 s
      1.427031      10303 <deep_copy>
      1.102441        531 cln_rhsaf
      1.067723       5000 vet_vgd_shell
      0.707400       7770 cln_pcr2
      0.596778        250 cln_drag
      0.488045        396 Kokkos::ViewCopy-4D
      0.471293        531 cln_fs0
      0.412539        339 cln_facecoef
      0.379074       1569 cln_apply_st
      0.316514         20 cln_gd_fs
      0.307091       2590 cls_restrict
      0.292541        636 cln_cslat
      0.260247        708 cln_fsch
      0.255582        137 cln_velocity
      0.254214       5080 vet_vgd_ring
```
