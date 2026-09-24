#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/vetcol2_0924/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/vetcol2_0924/gpu/log.err.%j
#SBATCH -J vetcol2_gpu
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# (1) GPU Cartesian bitwise gate of m1-vetcol2: ref = rt-integration 9c1a12d7 vs new.
# (2) vet_col build cost, team (default) vs the one-thread-per-column kernel, and the
#     total vs Eddington, He wedge grid 96 x 128 x 128, same binary, interleaved, 2 repeats.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/vetcol2_0924; G=$W/gpu; D=/viper/ptmp2/jinma/defaults_0923/cpu
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
A4=$W/inp2/sph_atm_vc.athinput
run a4t new none 1 -i $A4 mesh/nx1=64 meshblock/nx1=64 time/nlim=200 rad_m1/vet_col_surface_q=true
run a4s new none 1 -i $A4 mesh/nx1=64 meshblock/nx1=64 time/nlim=200 rad_m1/vet_col_surface_q=true rad_m1/vet_col_team=false
echo "== a4 team vs one-thread (GPU): $(python3 $W/cmp.py $G/a4t/new $G/a4s/new | tail -n1)"
python3 - $G/a4t/new $G/a4s/new <<'PY'
import glob, sys
sys.path.insert(0, '/viper/ptmp2/jinma/wt_vetcol2/vis/python')
import athena_read
fa = sorted(glob.glob(sys.argv[1] + '/tab/*.m1.*.tab'))[-1]
fb = fa.replace(sys.argv[1], sys.argv[2])
a = athena_read.tab(fa); b = athena_read.tab(fb)
import numpy as np
mx = 0.0
for k in a:
    if k.startswith('m1'):
        x, y = np.asarray(a[k]), np.asarray(b[k])
        mx = max(mx, float(np.max(np.abs(x - y)/np.maximum(np.abs(x), 1e-300))))
print('== a4 team vs one-thread max rel diff (last tab):', mx)
PY
H=$W/inp2/hewedge.athinput
for r in 1 2; do
  run hw_edd_$r new none 2 -i $H rad_m1/closure=eddington
  run hw_vct_$r new none 2 -i $H rad_m1/closure=vet_col
  run hw_vc1_$r new none 2 -i $H rad_m1/closure=vet_col rad_m1/vet_col_team=false
  run hw_vcq_$r new none 2 -i $H rad_m1/closure=vet_col rad_m1/vet_col_surface_q=true
done
for t in hw_edd_1 hw_vct_1 hw_vc1_1 hw_vcq_1 hw_edd_2 hw_vct_2 hw_vc1_2 hw_vcq_2; do
  echo "== $t: $(grep -h 'cpu time used\|zone-cycles/cpu_second\|vet_col: .*builds\|inner iterations' $G/$t/new/log.txt | tr '\n' ' ')"
done
echo GPU GATE DONE
