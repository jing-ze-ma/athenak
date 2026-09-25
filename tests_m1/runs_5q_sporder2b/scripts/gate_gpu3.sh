#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/sporder2b_0925/gpu/log3.out.%j
#SBATCH -e /viper/ptmp2/jinma/sporder2b_0925/gpu/log3.err.%j
#SBATCH -J sporder2b
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:10:00
# m1-sp-order2b GPU (3): bisect of the Cartesian GPU difference (snapshots new2, new3, new5)
# box2, boxh2: CPU bitwise, GPU not) rerun at TIGHT tolerance (implicit_tol 1e-12,
# implicit_lin_tol 1e-11): solver noise shrinks with the tolerance (gates.py rule 2),
# a logic change would not; plus ref vs ref (determinism)
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/sporder2b_0925; G=$W/gpu2; I=/viper/ptmp2/jinma/sph2_0924/inp
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
T="rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_lin_tol=1.0e-11"
run() {
  local t=$1 c=$2 p=$3 np=$4; shift 4
  local d=$G/$t/$c; rm -rf $d; mkdir -p $d; cd $d
  echo "== $t $c $(date +%T)"
  srun -n $np $W/bin/athena_${c}_${p}_gpu -d $d "$@" > log.txt 2>&1
  echo "rc=$?" >> log.txt
}
for c in new2 new3 new5; do
  run slab1t $c box 1 -i $I/slab2d_old_be.athinput time/tlim=100 problem/vpert=1.0e-2 $T
  run box2t $c box 2 -i $I/box3d_nd_be.athinput $B3 $T
done
echo GPU2 DONE
