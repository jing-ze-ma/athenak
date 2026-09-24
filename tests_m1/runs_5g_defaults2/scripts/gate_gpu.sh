#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/m1def2_0924/gpug/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/m1def2_0924/gpug/log.err.%j
#SBATCH -J m1def2g
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# m1-defaults2 GPU gate: be named (3 Cartesian cases) base 1560d07a vs new bitwise;
# key absent (new) vs hesdirk2 + implicit_vimp named (base) bitwise
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/m1def2_0924; G=$W/gpug; I=$W/inp
md5sum $W/bin/athena_base_box_gpu $W/bin/athena_new_box_gpu
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
run() {  # run <tag> <base|new> <np> <args...>
  local t=$1 c=$2 np=$3; shift 3
  local d=$G/$t/$c; rm -rf $d; mkdir -p $d; cd $d
  echo "== $t $c $(date +%T)"
  srun -n $np $W/bin/athena_${c}_box_gpu -d $d "$@" > log.txt 2>&1
  echo "rc=$?" >> log.txt
}
for c in base new; do
  run slab1 $c 1 -i $I/slab2d_old_be.athinput time/tlim=100 problem/vpert=1.0e-2
  run box2 $c 2 -i $I/box3d_nd_be.athinput $B3
  run boxvet2 $c 2 -i $I/box3d_nd_be.athinput $B3 rad_m1/closure=vet_sc time/nlim=30
done
run boxh2 new 2 -i $I/box3d_nd.athinput $B3 time/nlim=30
run boxh2 base 2 -i $I/box3d_nd_h2.athinput $B3 time/nlim=30
for t in slab1 box2 boxvet2 boxh2; do
  echo "== $t: $(python3 $W/scripts/gcmp.py $G/$t/base $G/$t/new | tail -n1) $(tail -n1 $G/$t/base/log.txt) $(tail -n1 $G/$t/new/log.txt)"
done
echo GPU GATE DONE
