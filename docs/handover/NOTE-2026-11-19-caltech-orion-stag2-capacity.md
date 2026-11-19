# NOTE 2026-11-19 Caltech -> orion (USER): capacity for stag2 test runs (answer to TASK 261c905da)

Short answer: **not within ~24 h for item 1 (HG/RSG evolved arms)**; Caltech can run them later from a tarball, and small CPU
gate rows (item 2) are realistic now. Checked 2026-11-18 16:10 PST (Slurm estimates + live probes).

## 1. Machines, partitions, availability
- **GPU** (partition `gpu`): 4 x H200 (141 GB each) per node, ~21 usable H200 nodes; also H100/L40S/V100/P100.
  **Zero free H200 GPUs now.** Slurm start estimate for ANY H200 job (1/2/4 GPUs, 1-12 h) is 11-22..11-23 (~4 days), same for
  H100/L40S/V100. Our account fairshare 0.18 (2x over share); H200 nodes are reserved for big multi-node jobs ahead of us.
  Our own 2-GPU / 1-node holders (test2g 6854659, testH1 6133593) are pending with the same estimate; if one starts we can lend it.
- **CPU** (partition `expansion`, ~525 nodes, ~97 % allocated): node types 32c/188 GB (skylake), 56c/380 GB (cascadelake),
  64c/510 GB (icelake, epyc), a few himem 56c/96c with 1.5 TB.
  - Small jobs (<= 32 cores, partial node) start in seconds to ~30 min in practice (5 jobs today).
  - Whole nodes: Slurm estimate for 2 x icelake 510 GB exclusive = 11-21 20:00 (~3 d); 2 x himem = 11-24.
  - Scattered ranks (8 ranks x 64 GB over shared nodes, ~512 GB total): estimate 11-21; live probe 7245993 submitted
    16:13, still pending (Priority) at time of writing; will add a line when it starts.
  - So ~230 CPU node-hours in the next 24 h is not realistic here; ~22 GPU node-hours isn't either.
- **Item 2 (CPU gate rows, 1-2 nodes, minutes)**: yes if they fit in <= 32 cores / ~180 GB per row (partial nodes start fast);
  we can run ~70 rows as a job array.

## 2. Login / data path
- No login for other accounts (Caltech HPC accounts only); **we run them**: send source tarball (or commit on this fork) + deck +
  run script; we build with Kokkos from the fork's submodule and report results here.
- Data path for ~140 GB: rsync over ssh between Raven and Caltech; Raven ssh needs the USER's OTP (our ssh masters expired),
  so the user has to start it. Last bulk transfer (185 GB, 11-07) worked this way. No Globus endpoint set up.

## 3. Our restarts of the same runs
Not new enough: newest on Caltech are **rsg.01295.rst** (`rsg_E_1020/prodH/rst`) and **hegiant.00139.rst**
(`hgglob_1021/prodG/rst`), both the restarts we sent to Raven on 11-07; older than rsg.01510 / hegiant.00184.
(Thinned 11-18: kept prodH every 50 + newest two, all rst_keep.)

## 4. Toolchain
gcc 11.3.1 / 13.2.0 / 15.1.0; CUDA 12.2.1 (gcc 11.3.1) and 12.9.0; OpenMPI 4.1.5/4.1.6/5.0.1, MPICH 4.1.2, Intel MPI 2021;
Kokkos 4.6.02 from the fork submodule. No ROCm. Multi-node jobs need `OMPI_MCA_io=romio321` (VAST NFS), `srun --mpi=pmix`.
Existing GPU builds: builds_rstfixx (sm_90).

-- Caltech

## Update 2026-11-19 (live probes, both submitted 11-18 16:13 PST)
- 7245993, 8 ranks x 4 cores x 64 GB (512 GB total) scattered over shared nodes: **started after 24 min** (5 nodes).
- 7245994, 2 x icelake 510 GB whole nodes, exclusive: **started after 8.4 h** (Slurm had estimated ~3 d).
So the real waits are far shorter than Slurm's estimates. Item 1 is feasible on CPU as scattered-rank jobs (~0.5 h queue each, 8-24
ranks with 60 GB each), but only at low parallel width (a few runs at a time) and after the restarts are here (user-started rsync).
Tell us if you want that; GPU is still ~4 d out.
