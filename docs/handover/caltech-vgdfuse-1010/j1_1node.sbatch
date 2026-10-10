#!/bin/bash
#SBATCH -J vgdfuse1
#SBATCH -A carnegie_poc
#SBATCH -p gpu
#SBATCH --exclude=hpc-sm-01-09
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:nvidia_h200:4
#SBATCH --cpus-per-task=8
#SBATCH --mem=480G
#SBATCH -t 01:00:00
#SBATCH -o /resnick/groups/carnegie_poc/jingze/vgdfuse_1010/j1.%j.out
# TASK-2026-10-10-caltech-vgdfuse-scaling, job 1 (1 node x 4 H200). Edit the paths in the CONFIG block only.
# Stops at the first arm with rc != 0 or FATAL (except the "may OOM" arms, which are reported and skipped).
module purge; module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
# ---- CONFIG ----
V=/resnick/groups/carnegie_poc/jingze/vgdfuse_1010
BB=/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_5b304cf9_nofma     # BASE
BN=/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_712ace27_nofma     # NEW
IA=/resnick/groups/carnegie_poc/jingze/agcarb_1009/files/agcar_rcxA_ge_accel_st_mg_hr.athinput
IB=/resnick/groups/carnegie_poc/jingze/bsg_1009/files/bsg_hr_dc5.athinput
export GPUWRAP=/resnick/groups/carnegie_poc/jingze/bsg_1009/gpuwrap.sh
export KTLIB=/resnick/groups/carnegie_poc/jingze/prof_1009/kt/profiling/simple-kernel-timer/kp_kernel_timer.so
R=$(dirname $(readlink -f $0))/vgdfuse_run.sh
[ -f "$R" ] || R=$V/vgdfuse_run.sh
# ----------------
S=$V/runs/$SLURM_JOB_ID; mkdir -p $S; nvidia-smi -L; hostname; md5sum $BB $BN $IA $IB
GR="mesh/nx2=128 mesh/nx3=128 mesh/x2min=1.3089969389957472 mesh/x2max=1.832595714594046 mesh/x3max=0.5235987755982988 meshblock/nx2=64 meshblock/nx3=64"
BHR="mesh/nx2=256 mesh/nx3=256 meshblock/nx2=64 meshblock/nx3=64"     # BSG hr: 16 blocks
go(){ # name bin input nranks nlim keys [mayoom]
  bash $R $S "$@" | tee -a $S/results.txt
  l=$(grep "^ARM $1 " $S/results.txt | tail -1)
  if ! echo "$l" | grep -q "rc=0 fatal=0"; then
    if [ -n "$7" ] && ! echo "$l" | grep -q "oom=0"; then echo "NOTE: $1 out of memory (allowed)"; return; fi
    echo "STOP: $1 failed"; exit 1
  fi; }
# (1) timing, 4 H200 = 1 MeshBlock per GPU, interleaved, 2 repeats: BASE, NEW (bitwise defaults),
#     NEW + lagged exchange (vet_gd_fuse_shells 2, 4; opt-in, NOT bitwise, failed the accuracy gate on A100)
for r in 1 2; do
  for c in agc bsg; do
    if [ $c = agc ]; then I=$IA; K=""; else I=$IB; K="$GR"; fi
    go ${c}_base_$r $BB $I 4 30 "$K"
    go ${c}_new_$r $BN $I 4 30 "$K"
    go ${c}_h2_$r $BN $I 4 30 "$K rad_m1/vet_gd_fuse_shells=2"
    go ${c}_h4_$r $BN $I 4 30 "$K rad_m1/vet_gd_fuse_shells=4"
  done
done
# kernel timer, one of each (KOKKOS_TOOLS_TIMER_BINARY=1: one .dat per rank; run kp_reader on them)
KT=1 go agc_base_kt $BB $IA 4 30 ""
KT=1 go agc_new_kt $BN $IA 4 30 ""
KT=1 go bsg_base_kt $BB $IB 4 30 "$GR"
KT=1 go bsg_new_kt $BN $IB 4 30 "$GR"
# (2) strong scaling AG Car A hr (4 blocks) on 1 and 2 H200 (4 = the timing arms above)
for n in 2 1; do
  go agc_s${n}_base $BB $IA $n 30 "" mayoom
  go agc_s${n}_new $BN $IA $n 30 "" mayoom
done
# (3) BSG hr full 16 blocks on 4 H200 (= strong point and the weak point at 4 blocks/GPU) and on 2 H200
go bsghr_s4_base $BB $IB 4 20 "$BHR" mayoom
go bsghr_s4_new $BN $IB 4 20 "$BHR" mayoom
go bsghr_s2_new $BN $IB 2 20 "$BHR" mayoom
go bsghr_s2_lm_new $BN $IB 2 20 "$BHR rad_m1/vet_gd_twin_lowmem=true" mayoom
# (4) weak scaling: BSG reduced (4 blocks) on 1 H200 = 4 blocks per GPU (vs bsghr_s4_* on 4 H200)
go bsg_w1_base $BB $IB 1 20 "$GR" mayoom
go bsg_w1_new $BN $IB 1 20 "$GR" mayoom
for f in $S/*_kt/*.dat; do [ -f "$f" ] && $(dirname $KTLIB)/kp_reader $f > ${f%.dat}.txt 2>&1; done
echo DONE $(date +%T)
