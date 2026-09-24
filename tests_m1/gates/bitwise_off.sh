#!/bin/bash -l
# The switch-off gate: implicit_op_check = 0 (named) and absent must be BITWISE the base
# binary.  usage: bitwise_off.sh <base athena> <new athena> <run dir> <input dir>
# (input dir as for gates.py, plus slab2d_nd.athinput of bench/m1_defaults_0923, which
# does not name the key).  Prints one cmp.py line per pair.
B=$1; N=$2; R=$3; I=$4; H=$(dirname $(readlink -f $0))
source /etc/profile.d/modules.sh 2>/dev/null; module purge >/dev/null 2>&1
module load gcc/14 openmpi/5.0 >/dev/null 2>&1
BOX="box3d_be_x mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=30"
LEV="rad_m1/implicit_vimp_fold=true rad_m1/implicit_fast_kernels=true rad_m1/implicit_one_pass=4 rad_m1/implicit_predictor_order=2"
OVL="rad_m1/implicit_halo_mpi=true rad_m1/implicit_halo_overlap=true rad_m1/implicit_halo_ovl_faces=true"
declare -A C
C[box1]="1|$BOX $LEV"
C[box2o]="2|$BOX $LEV $OVL"
C[box4o]="4|${BOX/meshblock\/nx2=16/meshblock/nx2=8} $LEV $OVL"
C[slab1]="1|slab2d_plm_vimp_be_x time/nlim=100 $LEV"
C[slab2o]="2|slab2d_plm_vimp_be_x meshblock/nx2=16 time/nlim=100 $LEV $OVL"
C[slabv2]="2|slab2d_plm_vimp_be_x meshblock/nx2=16 time/nlim=100 $LEV rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
C[slabnd2]="2|slab2d_nd meshblock/nx2=16 time/tlim=200"
for c in "${!C[@]}"; do
  IFS='|' read -r np args <<< "${C[$c]}"
  set -- $args; inp=$1; shift
  for w in base new; do
    exe=$B; [ $w = new ] && exe=$N
    d=$R/${c}_$w; rm -rf $d; mkdir -p $d
    (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none $exe \
       -i $I/$inp.athinput -d $d problem/vpert=1.0e-2 "$@" > log.txt 2>&1) &
  done
done
wait
for c in $(echo "${!C[@]}" | tr ' ' '\n' | sort); do
  echo "$c: $(python3 $H/cmp.py $R/${c}_base $R/${c}_new)"
done
