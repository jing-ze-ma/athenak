#!/bin/bash
# MR accuracy on the CPU He box (84x32x32, 4 blocks, vet_sc full, mg 3, hesdirk2):
# reference = k1 at cfl 0.075 with tight tolerances; arms at cfl 0.3: k1 (IMEX), k2, k4, k8;
# k1 at cfl 0.15 (IMEX time error scale).  All to the same tlim.  Usage: mr_box_cpu.sh EXE TAG TLIM
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
EXE=$1; TAG=$2; TL=$3
W=/viper/ptmp2/jinma/fast4_0925; R=$W/cpu/$TAG; mkdir -p $R
G="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/tlim=$TL time/nlim=-1 output1/dt=1.0e-30 problem/vpert=1.0e-2"
TI="rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_lin_tol=1.0e-11"
arm() { local n=$1; shift; mkdir -p $R/$n; cd $R/$n
  nice mpirun --oversubscribe --bind-to none -np 4 $EXE -i $W/inp/box.athinput -d . $G "$@" > log.txt 2>&1 & }
arm ref time/cfl_number=0.075 $TI
arm k1 time/cfl_number=0.3
arm k1h time/cfl_number=0.15
arm k2 time/cfl_number=0.3 rad_m1/implicit_mr_every=2
arm k4 time/cfl_number=0.3 rad_m1/implicit_mr_every=4
arm k8 time/cfl_number=0.3 rad_m1/implicit_mr_every=8
wait
echo MRBOX DONE
