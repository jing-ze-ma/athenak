#!/bin/bash
#SBATCH -J p16gate1
#SBATCH -p apudev
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --constraint=apu
#SBATCH --cpus-per-task=24
#SBATCH -t 00:15:00
#SBATCH -o /viper/ptmp2/jinma/ccpush16/gates1_mi300a.%j.out
# push-16 MI300A gates, job 1: arms A (push-16, new defaults forced off; must be bitwise to
# push 15) and B (push-16 defaults); 1-rank decks one lane per GPU, then 4-rank decks
# 2 ranks per GPU.  Reference = push-15 runs_B (push-15 binary, the same decks) + job 2.
P=/viper/ptmp2/jinma/ccpush16
NEW=$P/bin/athena_cc_hip_push16_n432
source /etc/profile.d/modules.sh; module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0
export KOKKOS_MAP_DEVICE_ID_BY=mpi_rank HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1 OMP_NUM_THREADS=1
. $P/one.sh
ONE=(g1_nltek1_x0.5 g1_hcyl_x0.5 g1_cb_x0.5 g1_shd_x0.25 g1_rwp g1_rwm g3_rhse_sp2 g3_rhse_sp3 g3_pulsew g3_spthin g4_rhse_cs6 g4_cs_thin1 g4_cs_feed1 g6_sp1x1wb_base g6_thin1_base g6_thin1l_base g1_rwm_x1 g1_shd_x0.25_x1 g4_cs_feed1_x1 g7_thin1_rot)
FOUR=(g2_nltek1_x0.5 g2_hcyl_x0.5 g2_cb_x0.5 g2_shd_x0.25 g2_rwp g2_rwm g3_rhse_sp3_mb4 g4_rhse_cs24)
echo "start $(date)"; md5sum $NEW
( for t in ${ONE[@]}; do one $NEW $P/pkg_A $P/runs_A $t 1 0; done ) &
( for t in ${ONE[@]}; do one $NEW $P/pkg_ref $P/runs_B $t 1 1; done ) &
wait
for t in ${FOUR[@]}; do for a in A B; do
  case $a in A) k=pkg_A;; B) k=pkg_ref;; esac
  one $NEW $P/$k $P/runs_$a $t 4 all; done; done
echo "ALLDONE $(date)"
