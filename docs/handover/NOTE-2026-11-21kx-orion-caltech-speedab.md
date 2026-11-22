# NOTE 2026-11-21kx orion -> caltech: SPEEDAB — old production vs new production speed + memory, HG and RSG (USER)

The USER asks how the speed of the new production compares with the old production. These are timing pairs on TEST H200 nodes only;
production keeps running.

## Binaries (build both with CUDA 12.9 + g++13, never 12.2; record md5)
- NEW `B_p1` = bf4f917d4 (branch caltech-prod-bf4f917d4). This is the build you already made for the successor holders.
- OLD HG  `B_o18` = 1d1972805 (branch origin/prod-1118). This is what the old viper HG production ran (athena_he_gpu72_x_1d197280_p1118).
- OLD RSG `B_o09` = 864d120d4 ("prod-1109r", on GitHub, e.g. origin/accbook-1113). This is what the old viper RSG production ran
  (athena_he_gpu72_q_864d120d_p1109, cmd.13255123).
- Use the same Kokkos/CMake recipe as B_p1, with the he_star problem.
- **Do you already have a prod-1118 or prod-1109 CUDA build? If so, name its path, sha and md5 and use it.**

## Restarts
Use your nearest old-prod restarts (Orion coord: fine for speed): HG `hgglob_1021/prodG/rst/hegiant.00139.rst`, RSG
`rsg_E_1020/prodH/rst/rsg.01295.rst`. If you have anything newer from the same lineage, use the newest. Name the file and its start cycle
C0. Each pair uses the SAME file for both arms.

## Keys
The key files are in docs/handover/speedab_1121/:
- {hg,rsg}_old.keys = the old viper production command lines;
- {hg,rsg}_new.keys = the current production .keys.

In those files, replace TABLEDIR with your copies of the tables, and check the md5s:

| table | md5 |
|---|---|
| ic_rsg_g25_deep.txt | 2d7de3c8 |
| rosseland_mesakap_4c0e8b64d26a.txt | fa0f4f42 |
| planck_ext2_x0.653848_z0.019757_d38e5526afef.txt | a9a8d6d3 |
| ic_hegiant_14p5M_mk.txt | f9678c2a |
| rosseland_mesakap_afe2dd8da2bd.txt | e2868b75 |
| planck_ext2_x0.0009496_z0.0215304_416a61196c64.txt | 23ee543f |

Append these keys to every arm, in this order:
`time/nlim=<C0+22> time/tlim=1e30 time/ncycle_out=1 job/final_restart=false output5/dt=-1 output2/last_time=1e30
output3/last_time=1e30 output4/last_time=1e30 output6/last_time=1e30 output5/last_time=1e30`
That gives 22 cycles with no restart dump and no data dumps. Delete any rst afterwards.

Env, all arms: KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 VGD_XFE=0 (plus your usual OMPI_MCA_io=romio321, srun --mpi=pmix).

## Layout
Both decks have 24 MeshBlocks. The old viper production peak was 160-178 GB per node for 6 blocks, so they do not fit on one node.
Every arm: **8 ranks x 3 blocks on 8 H200 = 2 nodes x 4 GPUs, 1 rank per GPU.**
Run the old and new arms of a pair on the SAME two nodes, one after the other, with no other job on those GPUs. Pair order:
  SA_HG_old (B_o18 + hg_old)  -> SA_HG_new (B_p1 + hg_new)
  SA_RSG_old (B_o09 + rsg_old) -> SA_RSG_new (B_p1 + rsg_new)
Run a memory sampler during each arm (`nvidia-smi --query-gpu=index,memory.used --format=csv,noheader -l 2` on both nodes) and report
the peak MiB per GPU.

## Report (V215/<arm>/res.txt + one table)
For each arm:
- steady s/cycle = (elapsed at cycle C0+21 - elapsed at C0+1) / 20, which excludes the first and the last cycle;
- the main-loop time;
- the dt range;
- peak GPU memory per GPU (max over the 8 GPUs);
- for B_p1, solves/stage: the mean and max of 'solves N' in the 'CLEAN cycle' lines;
- for the old binaries, their Newton/passes counters if printed;
- CW, FATAL, nan counts;
- the nodes used and the binary md5.

Priority: after the running test arms, before W3. Small: 4 arms x about 10 min, 2 nodes.

-- orion (for the USER)
