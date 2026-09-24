#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/h2fast_0924/gpug/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/h2fast_0924/gpug/log.err.%j
#SBATCH -J h2fgate
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# m1-h2fast GPU Cartesian be bitwise gate (the runs_5e_vetcol2 three cases): base 9c1a12d7 vs v4
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/h2fast_0924; G=$W/gpug; D=/viper/ptmp2/jinma/defaults_0923/cpu
md5sum $W/bin/athena_base_gpu $W/bin/athena_v4_gpu
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
run() {  # run <tag> <base|v4> <np> <args...>
  local t=$1 c=$2 np=$3; shift 3
  local d=$G/$t/$c; rm -rf $d; mkdir -p $d; cd $d
  echo "== $t $c $(date +%T)"
  srun -n $np $W/bin/athena_${c}_gpu -d $d "$@" > log.txt 2>&1
  echo "rc=$?" >> log.txt
}
for c in base v4; do
  run slab1 $c 1 -i $D/slab2d_old.athinput time/tlim=100 problem/vpert=1.0e-2
  run box2 $c 2 -i $D/box3d_nd.athinput $B3
  run boxvet2 $c 2 -i $D/box3d_nd.athinput $B3 rad_m1/closure=vet_sc time/nlim=30
done
for t in slab1 box2 boxvet2; do
  echo "== $t: $(python3 $W/gcmp.py $G/$t/base $G/$t/v4 | tail -n1) $(tail -n1 $G/$t/base/log.txt) $(tail -n1 $G/$t/v4/log.txt)"
done
echo GPU GATE DONE
