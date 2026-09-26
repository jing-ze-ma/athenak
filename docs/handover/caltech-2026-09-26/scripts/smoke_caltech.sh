#!/bin/bash
# smoke_caltech.sh <tag>: submit a 20-cycle WASP-121b 1x smoke test (base arm, production keys),
# one H200 (gpu) and one CPU node (8 MPI ranks), both on the debug QOS (30 min).
# Run dirs: /resnick/scratch/jingze/smoke_<tag>/{gpu,cpu}
T=${1:-c0926}
S=/resnick/home/jingze/ATHENAK/athenak/docs/handover/caltech-2026-09-26/scripts
B=/resnick/home/jingze/ATHENAK/builds
W=/resnick/scratch/jingze/smoke_$T
mkdir -p $W/gpu $W/cpu
cd $W/gpu && PFX=$W/gpu BIN=$B/athena_${T}_gpu NLIM=20 VERB=true WALL=00:20:00 \
  sbatch -q debug -t 00:25:00 -J smk_gpu $S/run_caltech.sub base
# CPU: same body, 8 ranks on a CPU node; the mesh has 6 x 2 x 2 = 24 meshblocks
cd $W/cpu && PFX=$W/cpu BIN=$B/athena_${T}_cpu NLIM=20 VERB=true WALL=00:20:00 \
  sbatch -q debug -t 00:25:00 -J smk_cpu -p expansion --gres=none --ntasks-per-node=8 --cpus-per-task=1 \
  --export=ALL,CPURUN=1 $S/run_caltech.sub base
