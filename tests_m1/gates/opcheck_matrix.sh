#!/bin/bash -l
# The implicit_op_check combination matrix (tests_m1/gates/README.md, check 1).
# usage: opcheck_matrix.sh <athena binary (CPU, MPI, PROBLEM=box_convection)> <run dir>
#        <input dir holding box3d_be_x.athinput and slab2d_plm_vimp_be_x.athinput>
# Every legal combination of implicit_vimp / implicit_vimp_fold / implicit_halo_mpi /
# implicit_halo_overlap / implicit_halo_ovl_faces / implicit_krylov_dev on 1, 2 and 4
# ranks, 3-D box and 2-D slab, one cycle, the check on the first two solves
# (implicit_op_check = -2: report only).  Prints one line per run and a summary.
# GEOS limits the geometries (default "box slab boxv slabv"; v = closure vet_sc).
exe=$1; R=$2; I=$3
source /etc/profile.d/modules.sh 2>/dev/null; module purge >/dev/null 2>&1
module load gcc/14 openmpi/5.0 >/dev/null 2>&1
mkdir -p $R
BOX="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8"
declare -A GEO
GEO[box1]="box3d_be_x $BOX meshblock/nx2=16 meshblock/nx3=16"
GEO[box2]="box3d_be_x $BOX meshblock/nx2=16 meshblock/nx3=16"
GEO[box4]="box3d_be_x $BOX meshblock/nx2=8 meshblock/nx3=16"
GEO[slab1]="slab2d_plm_vimp_be_x meshblock/nx2=32"
GEO[slab2]="slab2d_plm_vimp_be_x meshblock/nx2=16"
GEO[slab4]="slab2d_plm_vimp_be_x meshblock/nx2=8"
VET="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"   # non-zero edge coefficients
for n in 1 2 4; do GEO[boxv$n]="${GEO[box$n]} $VET"; GEO[slabv$n]="${GEO[slab$n]} $VET"; done
job() {   # name np geo overrides...
  local n=$1 np=$2 g=$3; shift 3
  local d=$R/$n; rm -rf $d; mkdir -p $d
  set -- ${GEO[$g]} "$@"
  local inp=$1; shift
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none $exe \
     -i $I/$inp.athinput -d $d problem/vpert=1.0e-2 time/nlim=1 \
     rad_m1/implicit_op_check=-2 "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt)
}
list=()
for geo in ${GEOS:-box slab boxv slabv}; do
  for np in 1 2 4; do
    for vimp in 0 1; do
      folds="0"; [ $vimp = 1 ] && folds="0 1"
      for fold in $folds; do
        V="rad_m1/implicit_vimp=$([ $vimp = 1 ] && echo true || echo false)"
        V="$V rad_m1/implicit_vimp_fold=$([ $fold = 1 ] && echo true || echo false)"
        for hm in 0 1; do
          H="rad_m1/implicit_halo_mpi=$([ $hm = 1 ] && echo true || echo false)"
          if [ $np = 1 ]; then
            for kd in 0 2; do
              list+=("${geo}_n${np}_v${vimp}f${fold}_h${hm}o0a0_k${kd}|$np|$geo$np|$V $H \
rad_m1/implicit_halo_overlap=false rad_m1/implicit_krylov_dev=$kd")
            done
          else
            ovls="0"; [ $hm = 1 ] && ovls="0 1 2"
            for o in $ovls; do
              O="rad_m1/implicit_halo_overlap=$([ $o -ge 1 ] && echo true || echo false)"
              O="$O rad_m1/implicit_halo_ovl_faces=$([ $o = 2 ] && echo true || echo false)"
              oo=$([ $o -ge 1 ] && echo 1 || echo 0); ff=$([ $o = 2 ] && echo 1 || echo 0)
              list+=("${geo}_n${np}_v${vimp}f${fold}_h${hm}o${oo}a${ff}_k0|$np|$geo$np|$V $H $O")
            done
          fi
        done
      done
    done
  done
done
for e in "${list[@]}"; do
  IFS='|' read -r n np g ov <<< "$e"
  job $n $np $g $ov &
  while [ $(jobs -rp | wc -l) -ge 24 ]; do wait -n; done
done
wait
np=0; nf=0
for e in "${list[@]}"; do
  IFS='|' read -r n rest <<< "$e"
  L=$R/$n/log.txt
  s=$(grep -c "^M1OPCHK solve [0-9]*: PASS" $L); f=$(grep -c "^M1OPCHK solve [0-9]*: FAIL" $L)
  nv=$(grep -o "([0-9]* variants" $L | head -1 | tr -d '(')
  rc=$(grep -o "rc=[0-9]*" $L | tail -1)
  worst=$(grep "^  var" $L | awk '{print $4}' | sort -g | tail -1)
  printf "%-32s solves pass %s fail %s  %-12s worst dy/row %-9s %s\n" $n "$s" "$f" \
    "$nv" "$worst" "$rc"
  if [ "$f" != 0 ] || [ "$s" = 0 ] || [ "$rc" != "rc=0" ]; then nf=$((nf+1)); fi
  np=$((np+1))
done
echo "MATRIX: $np runs, $nf with a FAIL, no check or a non-zero exit"
