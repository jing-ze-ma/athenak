#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/s1_0924/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/s1_0924/gpu/log.err.%j
#SBATCH -J s1_gpu2
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# GPU gate v2: ref 90b01f2f, rr = ref rerun, new ced90b9a, interleaved.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/s1_0924; G=$W/gpu2; D=/viper/ptmp2/jinma/defaults_0923/cpu
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
run() {  # run <tag> <ref|new> <np> <args...>
  local t=$1 c=$2 np=$3; shift 3
  local d=$G/$t/$c; rm -rf $d; mkdir -p $d; cd $d
  echo "== $t $c $(date +%T)"
  srun -n $np $W/bin/athena_${c}_box_convection_gpu -d $d "$@" > log.txt 2>&1
  echo "rc=$?" >> log.txt
}
for c in ref new rr; do
  run slab1 $c 1 -i $D/slab2d_old.athinput time/tlim=100 problem/vpert=1.0e-2
  run box2 $c 2 -i $D/box3d_nd.athinput $B3
  run boxvet2 $c 2 -i $D/box3d_nd.athinput $B3 rad_m1/closure=vet_sc time/nlim=30
done
for t in slab1 box2 boxvet2; do
  echo "== $t ref-vs-new: $(python3 $W/cmp.py $G/$t/ref $G/$t/new | tail -n1) $(tail -n1 $G/$t/ref/log.txt) $(tail -n1 $G/$t/new/log.txt)"
  echo "== $t ref-vs-ref: $(python3 $W/cmp.py $G/$t/ref $G/$t/rr | tail -n1) $(tail -n1 $G/$t/rr/log.txt)"
done
echo GPU GATE DONE
