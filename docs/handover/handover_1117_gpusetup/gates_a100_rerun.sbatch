#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --gres=gpu:a100:4
#SBATCH --time=01:00:00
#SBATCH -J cc_p15_gr
#SBATCH -o /raven/ptmp/jinma/ccpush15/gates_a100_rerun.%j.out
# rerun of the steps of 32417007 that failed on cpu bind (single-GPU steps on GPUs 1-3)
P=/raven/ptmp/jinma/ccpush15
REF=/raven/ptmp/jinma/ccgpu_1116/bin/athena_cc_a100_head62d7fixfin2_n432
NEW=$P/bin/athena_cc_a100_push15_n432
module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0
export KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1
one() {  # one <bin> <pkg> <runroot> <tag> <gpu>
  local B=$1 K=$2 R=$3 tag=$4 g=$5 d=$3/$4
  [ -d $d ] && mv $d $d.old.$(date +%m%d%H%M%S); mkdir -p $d; cd $d || return
  cp $K/decks/$tag.athinput in.athinput
  cmd="srun --overlap -N1 -n 1 -c 18 --cpu-bind=none -J p15_$tag env CUDA_VISIBLE_DEVICES=$g"
  echo "timeout 30m $cmd $B -i in.athinput -d ." > cmd.txt
  local t0=$(date +%s.%N)
  timeout 30m $cmd $B -i in.athinput -d . < /dev/null > run.log 2>&1; local rc=$?
  echo "wall $(awk -v a=$t0 -v b=$(date +%s.%N) 'BEGIN{printf "%.2f", b-a}')" >> run.log; echo "rc=$rc" >> run.log
  echo "$(date +%T) $(basename $R) $tag gpu=$g rc=$rc warn=$(grep -c 'CLEAN WARNING' run.log)"
}
ONE=(g1_nltek1_x0.5 g1_hcyl_x0.5 g1_cb_x0.5 g1_shd_x0.25 g1_rwp g1_rwm g3_rhse_sp2 g3_rhse_sp3 g3_pulsew g3_spthin g4_rhse_cs6 g4_cs_thin1 g4_cs_feed1)
G6=(g6_sp1x1wb_base g6_thin1_base g6_thin1l_base)
echo "start $(date) $(hostname)"; md5sum $REF $NEW
( for t in ${ONE[@]}; do one $NEW $P/pkg_A $P/a_runs_A $t 1; done ) &
( for t in ${ONE[@]}; do one $NEW $P/pkg_ref $P/a_runs_B $t 2; done ) &
( for t in ${G6[@]}; do one $NEW $P/pkg_ref $P/a_runs_B $t 3; done
  for r in 1 2; do one $NEW $P/pkg_ref $P/a_time_B_g3_r$r g5_sp1x1q_clean 3
    one $REF $P/pkg_ref $P/a_time_ref_g3_r$r g5_sp1x1q_clean 3; done ) &
wait
echo "ALLDONE $(date)"
