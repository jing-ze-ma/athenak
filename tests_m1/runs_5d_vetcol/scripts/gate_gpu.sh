#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/vetcol_0924/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/vetcol_0924/gpu/log.err.%j
#SBATCH -J vetcol_gpu
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# (1) GPU Cartesian bitwise gate of m1-vetcol: ref = m1-sp2 (s2_0924 new binary) vs new.
# (2) vet_col cost vs Eddington (and M1) on the He wedge grid at 96 x 128 x 128, same
#     binary, interleaved arms, 2 repeats.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/vetcol_0924; G=$W/gpu; D=/viper/ptmp2/jinma/defaults_0923/cpu
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
run() {  # run <tag> <ref|new> <prob> <np> <args...>
  local t=$1 c=$2 p=$3 np=$4; shift 4
  local d=$G/$t/$c; rm -rf $d; mkdir -p $d; cd $d
  echo "== $t $c $(date +%T)"
  srun -n $np $W/bin/athena_${c}_${p}_gpu -d $d "$@" > log.txt 2>&1
  echo "rc=$?" >> log.txt
}
for c in ref new; do
  run slab1 $c box_convection 1 -i $D/slab2d_old.athinput time/tlim=100 problem/vpert=1.0e-2
  run box2 $c box_convection 2 -i $D/box3d_nd.athinput $B3
  run boxvet2 $c box_convection 2 -i $D/box3d_nd.athinput $B3 rad_m1/closure=vet_sc time/nlim=30
done
for t in slab1 box2 boxvet2; do
  echo "== $t: $(python3 $W/cmp.py $G/$t/ref $G/$t/new | tail -n1) $(tail -n1 $G/$t/ref/log.txt) $(tail -n1 $G/$t/new/log.txt)"
done
H=$W/inp/hewedge.athinput
for r in 1 2; do
  run hw_edd_$r new none 2 -i $H rad_m1/closure=eddington
  run hw_vc_$r new none 2 -i $H rad_m1/closure=vet_col
  run hw_m1_$r new none 2 -i $H rad_m1/closure=m1
done
for t in hw_edd_1 hw_vc_1 hw_m1_1 hw_edd_2 hw_vc_2 hw_m1_2; do
  echo "== $t: $(grep -h 'cpu time used\|zone-cycles/cpu_second\|vet_col: .*builds\|inner iterations' $G/$t/new/log.txt | tr '\n' ' ')"
done
echo GPU GATE DONE
