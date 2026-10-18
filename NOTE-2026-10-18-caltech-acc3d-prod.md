# NOTE Caltech -> Raven (2026-10-18): 3-D accretor production copy (TASK-2026-10-18-caltech-acc3d-prod), status 1

All jobs SUBMITTED (Caltech, account carnegie_poc, --exclude=hpc-sm-01-09). Link 1 has NOT started yet; the START time
follows in the next NOTE update.

- Build: accretor-3d-1017 @ 186306cb, PROBLEM=ry_per_accretor, CUDA HOPPER90, MPI, Release, -ffp-contract=off (nofma),
  nvcc_wrapper, gcc 13.2 + cuda 12.9 + hpcx 2.17.1 (git-archive snapshot build, bsg_1009/build_inc_nofma.sh,
  BUILD_ROOT builds_accr). Two build jobs (4970290 expansion, 4970319 gpu, whichever runs first; lock-protected).
  Binary: /resnick/groups/carnegie_poc/jingze/builds_accr/bin/athena_gpu_ry_per_accretor_186306cb_nofma (md5 follows).
- Files: /resnick/groups/carnegie_poc/jingze/acc3d_prod_1018/files, all 5 md5s match the TASK (e9dab743, d81b7a46,
  99d668f6, 5b813886, 1ad77aed); deck acc3d_prod_mk.athinput (md5 9748aeaa) = .in with only the 3 @FILES@ paths substituted.
- Layout: 1 node x 4 H200, 4 ranks (srun --mpi=pmix -n 4 -c 8 --cpu-bind=cores gpuwrap.sh, CUDA_VISIBLE_DEVICES=SLURM_LOCALID,
  UCX_TLS=^xpmem), deck MeshBlocks unchanged (8 MBs, 2 per GPU).
- Pre-check 4970356 (A nlim 10 / B restart to 20 / C fresh 20, sect. 4) afterok on the build; it exits 1 unless
  B vs C bins differ in <= 1 byte, hydro.hst rows identical and user.hst cols 1-6,16-20,23 identical.
- Chain: link 1 4970357 afterok:4970356, link 2 4970358 afterok:4970357, both --kill-on-invalid-dep=yes,
  walltime 7-00:00:00 (normal QOS max; WT = 603600 s), acc3d_chain.sh unchanged; job rc = chain rc.
- Run path (rsync-able): /resnick/groups/carnegie_poc/jingze/acc3d_prod_1018 (static/, stream/, chain.txt, starts.txt).

## Update 2 (2026-10-18)
- Build 4970290 COMPLETED (5 min): binary md5 **9283c619823898f3d0196a3abca6129f**
  (/resnick/groups/carnegie_poc/jingze/builds_accr/bin/athena_gpu_ry_per_accretor_186306cb_nofma); 4970319 cancelled (redundant).
- Pre-check 4970356 TimeLimit cut to 40 min (backfill). Not started yet.
- USER GO: parallel **2-H200 variant** (DEVIATES from the TASK 4-GPU layout: 2 ranks, 4 MeshBlocks per GPU, same binary,
  same deck md5 9748aeaa) in /resnick/groups/carnegie_poc/jingze/acc3d_prod_1018/g2: pre-check 4972809 (40 min),
  link 1 4972811, link 2 4972812 (7 d, afterok + kill-on-invalid-dep). Whichever variant's link 1 STARTS first runs;
  the other set is cancelled at once. Start time + layout follow in the next update.
