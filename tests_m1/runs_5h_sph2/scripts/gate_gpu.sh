#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/sph2_0924/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/sph2_0924/gpu/log.err.%j
#SBATCH -J sph2gpu
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# m1-sph2 GPU: (1) Cartesian bitwise gate ref (d9583121 code) vs new: be named (3 cases)
# and key absent = the current hesdirk2 + vimp default (box3d_nd); (2) He wedge grid
# 96 x 128 x 128 (16 blocks on 2 GPUs, gas held, vet_col, 40 steps): be vs hesdirk2 per
# step, same binary, interleaved, 2 repeats
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/sph2_0924; G=$W/gpu; I=$W/inp
md5sum $W/bin/athena_ref_box_gpu $W/bin/athena_new_box_gpu $W/bin/athena_new_none_gpu
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
run() {  # run <tag> <ref|new> <box|none> <np> <args...>
  local t=$1 c=$2 p=$3 np=$4; shift 4
  local d=$G/$t/$c; rm -rf $d; mkdir -p $d; cd $d
  echo "== $t $c $(date +%T)"
  srun -n $np $W/bin/athena_${c}_${p}_gpu -d $d "$@" > log.txt 2>&1
  echo "rc=$?" >> log.txt
}
for c in ref new; do
  run slab1 $c box 1 -i $I/slab2d_old_be.athinput time/tlim=100 problem/vpert=1.0e-2
  run box2 $c box 2 -i $I/box3d_nd_be.athinput $B3
  run boxvet2 $c box 2 -i $I/box3d_nd_be.athinput $B3 rad_m1/closure=vet_sc time/nlim=30
  run boxh2 $c box 2 -i $I/box3d_nd.athinput $B3 time/nlim=30
done
for t in slab1 box2 boxvet2 boxh2; do
  echo "== $t: $(python3 $W/scripts/gcmp.py $G/$t/ref $G/$t/new | tail -n1) $(tail -n1 $G/$t/ref/log.txt) $(tail -n1 $G/$t/new/log.txt)"
done
for r in 1 2; do
  if [ $r = 1 ]; then A="be h2"; else A="h2 be"; fi
  for a in $A; do run hw_${a}_$r new none 2 -i $I/hwc_$a.athinput; done
done
for a in be h2; do for r in 1 2; do d=$G/hw_${a}_$r/new
  echo "== hw_${a}_$r: $(grep -h 'cpu time used\|zone-cycles/cpu_second' $d/log.txt | tr '\n' ' ') $(grep -h 'time_scheme=hesdirk2:\|vet_col: .*builds\|inner iterations\|NON-CONV' $d/log.txt | sort | uniq -c | tr '\n' ' ') $(tail -n1 $d/log.txt)"
done; done
echo GPU DONE
