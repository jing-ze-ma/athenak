#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --gres=gpu:a100:4
#SBATCH --time=02:00:00
#SBATCH -J cc_p16_gates
#SBATCH -o /raven/ptmp/jinma/ccpush16/gates_a100.%j.out
# push-16 A100 gates: build c249955c4, then ref = push-15 binary (default decks), A = push 16 with
# the new defaults forced off (bitwise to ref), B = push-16 defaults; bisect of the cs24 cycle-0
# CLEAN WARNING (one default off at a time); s/cycle + memory ref vs B interleaved
P=/raven/ptmp/jinma/ccpush16
REF=/raven/ptmp/jinma/ccpush15/bin/athena_cc_a100_push15_n432
NEW=$P/bin/athena_cc_a100_push16_n432
[ -e $NEW ] && mv $NEW $NEW.old.$(date +%m%d%H%M%S)
srun --overlap -N1 -n1 -c 72 --cpu-bind=none -J p16_build $P/build.sh || exit 1
[ -x $NEW ] || { echo NO_BINARY; exit 1; }
module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0
export KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1
one() {  # one <bin> <pkg> <runroot> <tag> <gpu 0-3|all>
  local B=$1 K=$2 R=$3 tag=$4 g=$5 d=$3/$4 q
  [ -d $d ] && mv $d $d.old.$(date +%m%d%H%M%S); mkdir -p $d; cd $d || return
  cp $K/decks/$tag.athinput in.athinput
  if [ $g = all ]; then
    cmd="srun --overlap -N1 -n 4 --ntasks-per-node=4 -c 18 --cpu-bind=none -J p16_$tag"; q=0,1,2,3
  else
    cmd="srun --overlap -N1 -n 1 -c 18 --cpu-bind=none -J p16_$tag env CUDA_VISIBLE_DEVICES=$g"; q=$g
  fi
  echo "timeout 20m $cmd $B -i in.athinput -d ." > cmd.txt
  ( while true; do echo "$(date +%T) $(nvidia-smi -i $q --query-gpu=index,memory.used --format=csv,noheader 2>/dev/null | tr '\n' ' ')" >> mem.log; sleep 1; done ) &
  local mp=$!
  local t0=$(date +%s.%N)
  timeout 20m $cmd $B -i in.athinput -d . < /dev/null > run.log 2>&1; local rc=$?
  kill $mp 2>/dev/null
  echo "wall $(awk -v a=$t0 -v b=$(date +%s.%N) 'BEGIN{printf "%.2f", b-a}')" >> run.log; echo "rc=$rc" >> run.log
  echo "$(date +%T) $(basename $R) $tag gpu=$g rc=$rc fatal/nan=$(grep -c -E '### FATAL|(^|[^A-Za-z_])-?(nan|NaN)([^A-Za-z_]|$)|DIVERGED' run.log) warn=$(grep -c 'CLEAN WARNING' run.log) retry=$(grep -c '^DT_RETRY' run.log)"
}
ONE=(g1_nltek1_x0.5 g1_hcyl_x0.5 g1_cb_x0.5 g1_shd_x0.25 g1_rwp g1_rwm g3_rhse_sp2 g3_rhse_sp3 g3_pulsew g3_spthin g4_rhse_cs6 g4_cs_thin1 g4_cs_feed1 g6_sp1x1wb_base g6_thin1_base g6_thin1l_base g1_rwm_x1 g1_shd_x0.25_x1 g4_cs_feed1_x1 g7_thin1_rot)
FOUR=(g2_nltek1_x0.5 g2_hcyl_x0.5 g2_cb_x0.5 g2_shd_x0.25 g2_rwp g2_rwm g3_rhse_sp3_mb4 g4_rhse_cs24)
T=(g5_sp1x1q_clean g4_cs_thin1 g4_cs_feed1_x1 g7_thin1_rot)
echo "start $(date) $(hostname)"; md5sum $REF $NEW
( for t in ${ONE[@]}; do one $REF $P/pkg_ref $P/runs_ref $t 0; done ) &
( for t in ${ONE[@]}; do one $NEW $P/pkg_A $P/runs_A $t 1; done ) &
( for t in ${ONE[@]}; do one $NEW $P/pkg_ref $P/runs_B $t 2; done ) &
( for k in ge yr fj; do one $NEW $P/pkg_bis $P/bis g4_rhse_cs6_$k 3; done ) &
wait
for t in ${FOUR[@]}; do one $REF $P/pkg_ref $P/runs_ref $t all; one $NEW $P/pkg_A $P/runs_A $t all
  one $NEW $P/pkg_ref $P/runs_B $t all; done
for k in ge yr fj; do one $NEW $P/pkg_bis $P/bis g4_rhse_cs24_$k all; done
for g in 0 1 2 3; do ( for t in ${T[@]}; do
  if [ $((g%2)) = 0 ]; then one $REF $P/pkg_ref $P/time_ref_g$g $t $g; one $NEW $P/pkg_ref $P/time_B_g$g $t $g
  else one $NEW $P/pkg_ref $P/time_B_g$g $t $g; one $REF $P/pkg_ref $P/time_ref_g$g $t $g; fi; done ) & done
wait
echo "ALLDONE $(date)"
