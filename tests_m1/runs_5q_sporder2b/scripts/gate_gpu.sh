#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/sporder2b_0925/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/sporder2b_0925/gpu/log.err.%j
#SBATCH -J sporder2b
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# m1-sp-order2b GPU: (1) bitwise ref (fe12a518) vs new: Cartesian be named (slab 1 rank,
# box3d_nd 2 ranks), key absent (box3d_nd, hesdirk2 + vimp), the He wedge with the OLD
# sp keys named (96 x 128 x 128, 16 blocks on 2 GPUs, 20 steps); (2) correctness: the
# T-S4 vet_col wedge (n = 64, 100 steps, hesdirk2, new sp defaults) and the same with a
# reflecting top, 1 GPU, compared with CPU afterwards; (3) cost on the He wedge grid:
# old keys vs new defaults with rbgs_fwd (pre) vs new defaults (mg_gc 1) vs
# time2_vet_col = rebuild, same binary,
# interleaved, 2 repeats, 40 steps
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/sporder2b_0925; G=$W/gpu; I=/viper/ptmp2/jinma/sph2_0924/inp; J=$W/inp
S=/viper/ptmp2/jinma/wt_sporder2b/tests_m1/runs_5q_sporder2b/scripts
md5sum $W/bin/athena_{ref,new}_{box,none}_gpu
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
run() {  # run <tag> <ref|new> <box|none> <np> <args...>
  local t=$1 c=$2 p=$3 np=$4; shift 4
  local d=$G/$t/$c; rm -rf $d; mkdir -p $d; cd $d
  echo "== $t $c $(date +%T)"
  srun -n $np $W/bin/athena_${c}_${p}_gpu -d $d "$@" > log.txt 2>&1
  echo "rc=$?" >> log.txt
}
run slab1 ref box 1 -i $I/slab2d_old_be.athinput time/tlim=100 problem/vpert=1.0e-2
run slab1 new box 1 -i $I/slab2d_old_be.athinput time/tlim=100 problem/vpert=1.0e-2
run box2 ref box 2 -i $I/box3d_nd_be.athinput $B3
run box2 new box 2 -i $I/box3d_nd_be.athinput $B3
run boxh2 ref box 2 -i $I/box3d_nd.athinput $B3 time/nlim=30
run boxh2 new box 2 -i $I/box3d_nd.athinput $B3 time/nlim=30
run hwo ref none 2 -i $J/hw_def.athinput time/nlim=20
run hwo new none 2 -i $J/hw_old.athinput time/nlim=20
for t in slab1 box2 boxh2 hwo; do
  echo "== $t: $(python3 $S/gcmp.py $G/$t/ref $G/$t/new | tail -n1) $(tail -n1 $G/$t/ref/log.txt) $(tail -n1 $G/$t/new/log.txt)"
done
run ts4g new none 1 -i $J/ts4_h2.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=100 output1/dcycle=100 output2/dcycle=100
run ts4r new none 1 -i $J/ts4_h2.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=100 output1/dcycle=100 output2/dcycle=100 rad_m1/implicit_bc_x1max=reflect
for r in 1 2; do
  if [ $r = 1 ]; then A="old pre def reb"; else A="reb def pre old"; fi
  for a in $A; do run hw_${a}_$r new none 2 -i $J/hw_$a.athinput; done
done
for a in old pre def reb; do for r in 1 2; do d=$G/hw_${a}_$r/new
  echo "== hw_${a}_$r: $(grep -h 'cpu time used' $d/log.txt | tr '\n' ' ') $(grep -ho 'NON-CONVERGED=[0-9.e+]*\|inner iterations mean=[0-9.e+]*\|[0-9.e+-]* ms per build\|time2_vet_col=[a-z]*' $d/log.txt | tr '\n' ' ') $(tail -n1 $d/log.txt)"
done; done
echo GPU DONE
