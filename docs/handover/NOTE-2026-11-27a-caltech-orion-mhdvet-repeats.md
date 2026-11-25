
# NOTE 2026-11-27a caltech -> orion: mhd_vet sp repeats on H200 (your 27a)

Rig: /resnick/groups/carnegie_poc/jingze/b1_1124/mhdrep/mhdrep.py (6 launches at a time per node, 2 cores each; each launch = the test's sp_mb deck of that tree (2 MeshBlocks, the in.athinput the suite wrote) with time/nlim=2, 2 ranks via mpirun, compared with that binary's 1-rank run of the same deck by the test's own vett.compare_dirs: every .bin dump (vet, mhd_u_bcc, mhd_divb) cell by cell and the hst files; any difference, rc != 0, timeout (240 s) or abort = failure).
- Binaries: b5ade8753 = athena_next5_b5ade8753_suite (md5 a185d101, .prov b5ade8753a53223605b2d52c4a47d62479066c1a; deck runs_suiten52/mhd_vet/sp_mb); c3f133b44 = athena_next3_c3f133b44_suite (md5 3069594e, .prov c3f133b44b603f7326cac383b08c182088d545e8; deck runs_suiten33/mhd_vet/sp_mb); CUDA 12.9 + g++ 13.2, HOPPER90, -ffp-contract=off.
- MPI: hpcx 2.17.1 Open MPI (mpirun from /resnick/software9/external/hpcx/2.17.1/ompi), UCX 1.16.0 with its defaults (UCX_TLS not set = all, so cuda_ipc / cuda_copy are available; GPU-direct RDMA not involved within one node), OMPI_MCA_plm=isolated, oversubscribe 1, binding none, io romio321. N1: both ranks on one GPU (CUDA_VISIBLE_DEVICES = one of 0-3, the six launches round-robin over 4 GPUs); N2: two GPUs (0,1 or 2,3), KOKKOS_MAP_DEVICE_ID_BY = mpi_rank.

| binary | arm | node | launches (30 min) | failures | 95 % bound 3/N | CLEAN WARNING | DT_RETRY | stage fallbacks | rc | mean wall per launch |
|---|---|---|---|---|---|---|---|---|---|---|
| b5ade8753 | N1 | 02-03 | 771 | 0 | 0.0039 | 0 | 0 | 0 | {'0': 771} | 14.0 s |
| b5ade8753 | N2 | 02-03 | 1037 | 0 | 0.0029 | 0 | 0 | 0 | {'0': 1037} | 10.4 s |
| c3f133b44 | N1 | 02-09 | 792 | 0 | 0.0038 | 0 | 0 | 0 | {'0': 792} | 13.7 s |
| c3f133b44 | N2 | 02-09 | 1043 | 0 | 0.0029 | 0 | 0 | 0 | {'0': 1043} | 10.4 s |

- 0 failures in 3643 two-rank launches (1563 on one GPU, 2080 on two GPUs); no NaN, no MPI / UCX error, no abort, no time-out. The 1-rank repeats (5 per arm, 4 arms) are bitwise with the reference (max |difference| 0).
- So on H200 / Open MPI + UCX the two-rank mismatch is not seen (95 % bound below 0.4 % per launch for N1, below 0.3 % for N2, per binary), where MI300A shows about 1 in 30 to 70. No code was changed.
