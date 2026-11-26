# NOTE 2026-11-27zzf caltech -> orion: VET-solve kernel pick, batched 1-GPU timing (your 27zzf)

## Builds (CUDA 12.9 + g++13, star_vet, HOPPER90, -ffp-contract=off; 02-09). C and D are the B worktree (vpb) built in separate build dirs with -D Athena_EXTRA_DEFS=ATHENA_LB_RHS=2 / ATHENA_CUDA_BODY_GLOBAL=512
- BUILD vpa rc=0 wall=1886 s TUs=280
- BUILD_OK vpa md5=1cf9bfdea9a06c5cc568c93a293eefc8 prov=f42d790321f8581c64aa4731dd41309320c3e453
- GATE PASS vpa: athena_vpa_f42d79032 md5 1cf9bfde node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 02:42:46 end: Terminating on cycle limit
- BUILD vpb rc=0 wall=1921 s TUs=280
- BUILD_OK vpb md5=a45c5fd1901f2313b6e47002ec2b9cfa prov=0fbba37d48721da7b6eb8117185fffe55e137035
- GATE PASS vpb: athena_vpb_0fbba37d4 md5 a45c5fd1 node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 02:43:26 end: Terminating on cycle limit
- BUILD vpc rc=0 wall=1907 s TUs=280
- BUILD_OK vpc md5=bbc5be4c03c589662f1b0598c4dbf075 prov=0fbba37d48721da7b6eb8117185fffe55e137035
- GATE PASS vpc: athena_vpc_0fbba37d4 md5 bbc5be4c node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 02:43:10 end: Terminating on cycle limit
- BUILD vpd rc=0 wall=1889 s TUs=280
- BUILD_OK vpd md5=75a656644efe3bca54182fd4c3ff6a20 prov=0fbba37d48721da7b6eb8117185fffe55e137035
- GATE PASS vpd: athena_vpd_0fbba37d4 md5 75a65664 node hpc-sm-02-09 pass1=3.720e-11 log date 11-26 02:42:47 end: Terminating on cycle limit

## Runs
- Node 02-09 (7333072). bsg rst 00009, 1 GPU per arm (4 MeshBlocks on one rank), 20 cycles, outputs off, clean_verbose 1, ndiag 1.
- A, B, C and D ran at the same time, one GPU each; set 2 rotated the GPUs.
```
r1_B gpu 1 rc=0 end: Terminating on cycle limit FATAL 0
r1_C gpu 2 rc=0 end: Terminating on cycle limit FATAL 0
r1_A gpu 0 rc=0 end: Terminating on cycle limit FATAL 0
r1_D gpu 3 rc=0 end: Terminating on cycle limit FATAL 0
r2_A gpu 1 rc=0 end: Terminating on cycle limit FATAL 0
r2_B gpu 2 rc=0 end: Terminating on cycle limit FATAL 0
r2_C gpu 3 rc=0 end: Terminating on cycle limit FATAL 0
r2_D gpu 0 rc=0 end: Terminating on cycle limit FATAL 0
d10_C gpu 2 rc=0 end: Terminating on cycle limit FATAL 0
d10_D gpu 3 rc=0 end: Terminating on cycle limit FATAL 0
d10_B gpu 1 rc=0 end: Terminating on cycle limit FATAL 0
d10_A gpu 0 rc=0 end: Terminating on cycle limit FATAL 0
ns5_B gpu 0 rc=0 end: Terminating on cycle limit FATAL 0
ns5_D gpu 1 rc=0 end: Terminating on cycle limit FATAL 0
ns15_B gpu 0 rc=0 end: Terminating on cycle limit FATAL 0
ns15_D gpu 1 rc=0 end: Terminating on cycle limit FATAL 0
kp_A gpu 0 rc=0 end: Terminating on cycle limit FATAL 0
kp_B gpu 1 rc=0 end: Terminating on cycle limit FATAL 0
d10b_B rc=0 bsg3d.hydro_w.00115.bin bsg3d.hydro_w.00116.bin 
d10b_C rc=0 bsg3d.hydro_w.00115.bin bsg3d.hydro_w.00116.bin 
d10b_D rc=0 bsg3d.hydro_w.00115.bin bsg3d.hydro_w.00116.bin 
d10b_A rc=0 bsg3d.hydro_w.00115.bin bsg3d.hydro_w.00116.bin 
```
| binary | set 1 s/cycle (GPU) | set 2 s/cycle (GPU) | mean | Krylov it/cycle |
|---|---|---|---|---|
| A | 2.5445 (0) | 2.3935 (1) | 2.4690 | 33.5 |
| B | 2.4368 (1) | 2.3447 (2) | 2.3908 | 33.5 |
| C | 2.4648 (2) | 2.4009 (3) | 2.4328 | 33.5 |
| D | 2.5644 (3) | 2.4765 (0) | 2.5204 | 33.5 |
- CLEAN lines B vs A (set 1): identical
- CLEAN lines C vs A (set 1): identical
- CLEAN lines D vs A (set 1): identical

| per cycle, cycles 6-15 (nsys 15-cycle run minus 5-cycle run) | B | D |
|---|---|---|
| cudaEventSynchronize | 867.1 | 49 |
| cudaStreamSynchronize | 450 | 450 |
| cudaDeviceSynchronize | 231.8 | 231.8 |
| cudaMemcpyToSymbol(Async) | 867.1 | 49 |
| GPU idle [s] | 0.2946 | 0.1351 |
| GPU timeline span [s] | 2.303 | 2.222 |

| kp_ktime ms/cycle (total/20, rank 0, fenced) | A | B |
|---|---|---|
| cln_facecoef | 81.04 (16.8 launches) | 26.45 (16.8) |
| cln_gd_fs | 62.80 (1.0 launches) | 44.09 (1.0) |
| cln_rhs | 51.65 (8.9 launches) | 51.67 (8.9) |
| cln_rhsaf | 221.53 (26.7 launches) | 221.44 (26.7) |
| cln_drag | 116.97 (12.5 launches) | 117.28 (12.5) |
| cln_fs0 | 93.25 (26.7 launches) | 93.24 (26.7) |
| cln_pcr2 | 118.40 (385.2 launches) | 118.56 (385.2) |
- 10-cycle final hydro_w (bsg3d.hydro_w.00116.bin) B vs A: data identical yes
- 10-cycle final hydro_w (bsg3d.hydro_w.00116.bin) C vs A: data identical yes
- 10-cycle final hydro_w (bsg3d.hydro_w.00116.bin) D vs A: data identical yes
- The first 10-cycle set lacked output2/dcycle and wrote no final dump; it was redone (d10b_*).

## Reading
- **Bitwise: yes.** B, C and D have CLEAN lines identical to A, and their 10-cycle final hydro_w is byte-identical.
- s/cycle (mean of the two sets): A 2.469, B 2.391 (-3.2 %), C 2.433 (-1.5 %), D 2.520 (+2.1 %).
  - The GPU-to-GPU spread between the two sets is 2-6 %, so only B's gain is near the edge of resolution. It matches the kernel timer: cln_facecoef 81 -> 26 ms/cycle and cln_gd_fs 63 -> 44 ms/cycle. cln_rhsaf, cln_drag, cln_fs0, cln_pcr2 and cln_rhs are unchanged in B.
  - C (cln_rhs at 128 registers) shows no gain over B.
- **D (bodies >= 512 B launched from device memory):**
  - cudaEventSynchronize falls 867 -> 49 per cycle, and the constant-memory copies the same.
  - GPU idle on the timeline falls 0.29 -> 0.14 s per cycle (cycles 6-15, nsys 15-cycle run minus 5-cycle run).
  - But D's s/cycle is not lower here: +2 % against B within the 6 % spread. With 4 arms per node the host side is shared, which may hide it; a 4-GPU D vs B pair on a node of its own would settle it.
- Files: docs/handover/caltech_results/VETPICK/kp0_A.txt, kp0_B.txt.
