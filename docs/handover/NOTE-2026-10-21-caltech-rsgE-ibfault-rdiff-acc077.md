# NOTE 2026-10-21 caltech -> raven: RSG E link 1 died (IB fault) + requeued; rdiff 563 result R != C; accretor 077dda23 build queued

## RSG E (stays on 4de8168d 9a56b4a1 + clip_rst5d a7e734c4, run.cfg unchanged)
- Link 1 5126654 started 10-21 01:47 PDT on hpc-sm-02-[04,11], died after 63 s, still inside ic_rad (IcRadRun -> ImplicitSolve ->
  VetColBuild -> VetGdBuild -> VetGdCsSource -> VetGdCsHalo -> MPI_Waitall): UCX "Transport retry count exceeded on mlx5_3:1/IB",
  RDMA_READ 5.78 MB from GPU memory, rank 0. No code error, no rst written. Same node pair ran the 2-node g25 det check 5008036 fine,
  but that deck had no cubed-sphere vet_gd halo, so this is the first 2-node run of the vet_gd CS band exchange here.
- Requeued with no changes: link 2 5126655 had its dependency removed and is now link 1 (pending, queue age kept) -> 5126656 afterok.
  Old log kept as prod/run.fail_5126654_ib.log. If it fails the same way again, next try = UCX_IB_GPU_DIRECT_RDMA=no (stage through host).
- DeltaAI's copy is the one actually running.

## rdiff 5140304 (1 x 4 H200, 9a56b4a1): NOT bitwise
- A (deck E N=32): R vs C DIFF in hydro.hst (first t=1.02e5), user.hst (t=1.01e5), hydro_w 00002/00003, m1 00002/00003, m1_vet 00001.
  R2 vs R DIFF too (first t=1.60e5 / 1.57e5): two restarts from the same rst disagree -> run-to-run nondeterminism, not only restart state.
- B (old rsg_det N=48 + vet_rotate_avg=1): RB vs CB DIFF (t=3.35e4), RB2 vs RB DIFF (t=5.22e4); all files DIFF.
- All 6 runs rc=0, FATAL 0. Dirs kept: /resnick/groups/carnegie_poc/jingze/rsg_E_1020/rdiff563/{C,R,R2,CB,RB,RB2}, log rd.5140304.out.

## Accretor (V12_READY)
- Build 077dda23 ry_per_accretor H200 nofma = job 5223087 (debug queue, pending). md5 goes in the next NOTE.
- Fresh run dir /resnick/groups/carnegie_poc/jingze/acc3d_v12_1021: deck from template 5de0635e (@FILES@ -> run dir files),
  acc3d_prod_v12.athinput md5 6b771411; check deck = the same + hst dcycle=1. Once the build is done: link1 (in-job C/C2/R check, then chain)
  + 5 links afterok, 2 H200, 7 d each.
