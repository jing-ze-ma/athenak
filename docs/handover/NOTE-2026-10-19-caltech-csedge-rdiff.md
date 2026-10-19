# NOTE 2026-10-19 caltech -> raven: CSEDGE arms + RSG restart-diff (R2/RV/CV) SUBMITTED (no results yet)

Answers TASK-2026-10-19-caltech-csedge and TASK-2026-10-19-caltech-rsg-restart-diff. Binary a2713be0 H200 build,
md5 e31ea39ea5e05148e2dbea0b95328b27 (builds_fin51/bin/athena_gpu_he_star_m1_a2713be0_nofma); deck = rsg_g25_prod.athinput
(b49f142d + local IC/table paths only).

## A) CSEDGE: one Slurm job per arm, 1 H200 each, 24 blocks (N=48 mb 24x24; B24: N=24 mb 24x24 = 6 blocks)
| arm | job | extra key |
|---|---|---|
| B | 5059457 | - |
| F | 5059458 | problem/mlt_closure=frozen |
| H | 5059459 | hydro/rsolver=hlle |
| N | 5059460 | rad_m1/gas_feedback=false |
| V | 5059461 | rad_m1/vet_gd=false |
| B24 | 5059462 | mesh/nx2=24 mesh/nx3=24 |

Command (script caltech/arm.sh; exact line saved per arm as csedge_1019/<arm>/cmdline.txt):
`srun --mpi=pmix -N 1 -n 1 -c 8 --cpu-bind=cores gpuwrap.sh $BIN -i rsg_g25_prod.athinput mesh/nx2=48 mesh/nx3=48 meshblock/nx2=24
meshblock/nx3=24 time/tlim=1.04e6 time/ndiag=100 output1/dt=2000 output2/dt=43200 output3/dt=43200 output4/dt=1e30 output5/dt=1e30
output6/dt=1e30 <arm key> -d csedge_1019/<arm>`
Walltime estimate: N=48 on 4 H200 ran ~0.6-0.8 s/cycle (det1, dt 329 s) -> 1.04e6 s ~ 3200 cycles; 1 H200 guess 2-3 s/cycle
-> ~2-3 h if dt holds (limit 16 h). Dirs: /resnick/groups/carnegie_poc/jingze/rsg_g25_prod/csedge_1019/<arm>/.

## B) restart diff: job 5059481 (1 node x 4 H200, 6 blocks/GPU, same layout as det1 job 5008037)
R2 = identical restart from n48/C rsg.00001.rst (cycle 100) to cycle 200; CV = continuous t=0 -> 200 with rad_m1/vet_rotate_avg=1;
RV = restart from the same rst with rad_m1/vet_rotate_avg=1. Compare (caltech/rdiff/compare.sh): R2 vs R, RV vs CV (+ R2 vs C,
CV vs C), hst rows + bin data after <par_end>. Dir /resnick/groups/carnegie_poc/jingze/rsg_g25_prod/det/rdiff_1019/.

Exact C and R command lines of the earlier run (det1.sh, job 5008037; copied as caltech/rdiff/det1.sh), deck det/rsg_det.athinput =
rsg_g25_prod.athinput + `dcycle` per output (output1 hst dcycle=1, output2 hydro_w 100, output3 m1 100, output4 log 50,
output5 rst 100; output6 m1_vet dcycle=100 but disabled by dt=1e30 on the command line):
- C: `srun --mpi=pmix -N 1 -n 4 --ntasks-per-node=4 -c 8 --cpu-bind=cores gpuwrap.sh $BIN -i det/rsg_det.athinput mesh/nx2=48 mesh/nx3=48 meshblock/nx2=24 meshblock/nx3=24 time/nlim=200 output6/dt=1e30 -d det/n48/C`
- R: `srun ... $BIN -r det/n48/C/rst/rsg.00001.rst mesh/nx2=48 mesh/nx3=48 meshblock/nx2=24 meshblock/nx3=24 time/nlim=200 output6/dt=1e30 -d det/n48/R`
(gpuwrap.sh: modules + CUDA_VISIBLE_DEVICES=$SLURM_LOCALID, UCX_TLS=^xpmem; OMPI_MCA_io=romio321 exported.)
Earlier result (det1): C == C2 bitwise; R vs C hst first diff t=4.633e4 (user.hst, cycle 141) / 4.764e4 (hydro.hst, cycle 145).

## Collect
`/resnick/groups/carnegie_poc/jingze/rsg_g25_prod/csedge_1019/collect.sh <bsg-files-1009 checkout>` runs csedge_diag.py on the dumps
so far and copies hst / run.log tail30 / cmdline / s_per_cycle / diag.txt to docs/handover/csedge-1019/caltech/<arm>/ and the
restart-diff compare to caltech/rdiff/compare.txt; then commit + push. Results follow in a later commit.
