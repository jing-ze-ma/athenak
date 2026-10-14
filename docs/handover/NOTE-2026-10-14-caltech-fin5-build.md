# NOTE caltech -> raven: FIN5 build + smoke on H200 (TASK-2026-10-14-caltech-fin5-build, a1f3a175)

Caltech, 2026-10-14. Work dir /resnick/groups/carnegie_poc/jingze/fin5_1014 (smoke.sh, build.sh, smk_*.out).
No running or queued production was touched.

## 1. Builds (compute jobs, expansion debug; same incremental tree, toolchain and flags as the FIN4 d505d21f / 8e018272 / 5b7e1bc5 nofma builds)
gcc 13.2.0 + cuda 12.9.0 + hpcx 2.17.1 ompi, Kokkos 4.6.2 (d8e9af03), `-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON
-D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS=-ffp-contract=off`,
default Athena_M1_VGD_NMAX, no -DM1_VGD_NCCL (no ACC_NCCL string in the binary).

| binary | git rev-parse HEAD | md5 | job |
|---|---|---|---|
| athena_gpu_he_star_m1_ff03587a_nofma | ff03587a8c11dedb8b00fcd72fed2fcd4cc03581 | c05170d98488f1f8fc4978c00b763fbc | 4527765 |
| athena_gpu_he_star_m1_7299d0d7_nofma | 7299d0d77d17cbc21a9ce378c74f67dcd25e17db | b71887f49368b556cc816761780c184b | 4527765 |
| athena_gpu_he_star_m1_300b6ce8_nofma (extra Linhao base) | 300b6ce82a2a8157dcbf0d2e10f3be62b8280747 | fb938d1eb44ac16ac14e11a4b5a436f9 | 4527864 |
| base HG: athena_gpu_he_star_m1_d505d21f_nofma | d505d21f | 1440c05320a7bef88c4ad960019b9840 | (FIN4) |
| base LI: athena_gpu_he_star_m1_8e018272_nofma | 8e018272 | 4a92728d016609effe6d4586fb2e55bf | (FIN4) |

Binaries in /resnick/home/jingze/ATHENAK/builds/bin.

## 2. Smoke (gpu debug QoS, 1 node x 2 H200, 2 ranks, production launch: srun --mpi=pmix -n 2 -c 8, 1 GPU per rank, UCX_TLS=^xpmem)
Jobs: He giant 4531633 (A B C) + 4531634 (D E); Linhao 4531635 (A B C) + 4531636 (D E F).
Decks: hegiant_fin4.athinput (f266cd9b), bsgLI_bc_caltech.athinput (d77f9a10). Keys: time/nlim (12; C = 6), time/ndiag=1,
output1 (hst) dt=1e-6 (every cycle), all other outputs dt=1e30 (written at t=0 and at the end, rst at the end).
Arms: A FIN4 base fresh 12 | B FIN5 fresh 12 | C FIN5 fresh 6 | D FIN5 `-r` C -> 12 | E FIN5 + VGD_RING_POISON=1 fresh 12 |
F (Linhao only) 300b6ce8 fresh 12.
rst compare = bytes after `<par_end>`; bin compare = bytes after `<par_end>` of every hydro_w / m1 / m1_vet (/ m1_face) dump;
hst = every data row.

| pair | He giant | Linhao bc |
|---|---|---|
| FIN5 vs FIN5-poison (B/E) | BITWISE: rst, all bins, user+hydro hst (13 rows) | BITWISE: rst, all bins, user+hydro hst (13 rows) |
| FIN5 vs FIN5 restart 6 -> 12 (B/D) | BITWISE: rst, final bins; all D (and C) hst rows verbatim in B | BITWISE: rst, final bins; all D (and C) hst rows verbatim in B |
| 300b6ce8 vs FIN5 (F/B) | - | BITWISE: rst byte-equal, all bins, both hst |
| FIN4 vs FIN5 (A/B) | DIFFERS from t = 0 (expected, see 2b) | data BITWISE; differs ONLY in the L_in column and the pgen trailer (expected, 300b6ce8) |

Linhao A/B detail: hydro.hst identical; user.hst identical except column 7 `L_in` (rows 1..12; row 0 identical);
all hydro_w / m1 / m1_vet / m1_face bins (t=0 and cycle 12) byte-equal; rst: the last 9 043 274 496 bytes (the whole
MeshBlock data block) byte-equal, the files differ only inside the first 676 bytes after `<par_end>` (restart headers)
and B is 16 bytes longer (the new L_in pgen trailer: Real L + tag). F (300b6ce8) vs B is fully byte-equal.

He giant A/B: rst data differ from byte 669 after `<par_end>`; hydro_w, m1, m1_vet bins differ already at t = 0; both hst differ.
The FIN4 A log has `ic_rad: no fused twin intensities`, FIN5 (B..E) do not -> the expected FIN4 IC overflow (sect. 2b),
not a FIN5 problem. FIN5 vs poison and FIN5 restart are bitwise.

Log checks, all arms: FATAL 0, NaN 0, NON-CONVERGED>0 0; `out-of-slab accesses main 0 twin 0` in every arm;
FIN5 arms `vet_gd_band_ring ... dense band saves 0` (E prints "(POISON test)").
Ring sizes rank 0: HG compact intensities 7.10 GB + 3 band slots 0.21 GB; LI 8.24 GB + 0.57 GB.

### s/cycle (cycles 3..11, one sample each) and GPU peak (nvidia-smi max, both GPUs equal)
| | FIN4 s/cyc | FIN5 s/cyc | FIN5-poison | FIN4 peak | FIN5 peak | memory |
|---|---|---|---|---|---|---|
| He giant | 0.954 (d505d21f) | 0.865 (-9.3 %) | 0.915 | 85 653 MiB (83.6 GiB) | 49 949 MiB (48.8 GiB) | -41.7 % |
| Linhao bc | 1.280 (8e018272); 1.334 (300b6ce8) | 1.254 (-2.1 % / -6.0 %) | 1.309 | 85 383 MiB (83.4 GiB) | 51 143 MiB (49.9 GiB) | -40.1 % |

Single 8-cycle samples on one node (hpc-sm-02-10), all arms sequential in the same jobs: speed differences are within the
expected few %; the memory drop is the clear result.

## 2b. FIN4 fresh-start logs: `ic_rad: no fused twin intensities` IS PRESENT
Line: `<rad_m1> ic_rad: no fused twin intensities: J is the raw formal solution (ray noise not removed)`, once, in
- /resnick/groups/carnegie_poc/jingze/hegiant_1014/fin4/run.log (He giant FIN4 production fresh start, chain hgF4),
- hegiant_1014/val/A/run.log and val/C/run.log (FIN4 validation),
- hegiant_large_1014/gate/A/run.log and gate/C/run.log (large-domain gate).
So every Caltech He giant FIN4 fresh start (normal and large domain) began with the degraded radiation IC (no twin J
correction). Note: the large-domain gate C, whose background run failed the E0 photosphere-band comparison
(NOTE 75437e36), also started from this degraded IC, so the IC skip may bear on that mismatch. Raven decides.
Nothing was restarted or changed.

## Not done / open
- FIN5 not used for any production; E1 stays on 5b7e1bc5 (FIN5 has no dt_retry).
- Run dirs hold ~10 GB rst per arm (fin5_1014/hg, fin5_1014/li); can be pruned once Raven has read this.
