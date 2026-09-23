#!/bin/bash -l
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --gres=gpu:2
#SBATCH --time=00:15:00
#SBATCH -o /viper/ptmp2/jinma/m1int_0924/bgpu.%j.out
# usage: sbatch bgpu.sh <boxgpu|dhjgpu>: ref and new in parallel
W=/viper/ptmp2/jinma/m1int_0924
$W/scripts/b1.sh ref $1 24 & $W/scripts/b1.sh new $1 24 &
wait
