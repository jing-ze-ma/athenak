#!/bin/bash -l
#SBATCH -o log.out.%j
#SBATCH -e log.err.%j
#SBATCH -J cksw_a1
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# A1: ck_spherical=false (plane-parallel), SAME binary as bench/prof_0922/plain.
# (rot 283, ncycle 4429831), production input (bench/ck_sph_ab/sph/deep_hot_jupiter.athinput
# with rt_cell_report/rt_report_every removed), 2 ranks/2 GPUs, ~300 cycles via nlim.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"
N=4430131   # ncycle(4429831) + 300
/usr/bin/time -f "PLAIN_WALL %e s" srun -n 2 ./athena \
     -r /viper/u2/jinma/ATHENAK/bench/cs_mhd_prod3/rst/dhj.00567.rst \
     -i deep_hot_jupiter.athinput \
     -t 00:13:00 \
     time/nlim=$N \
     > plain.log 2>plain.err
echo "rc=$?"
grep -E "PLAIN_WALL|cpu time used|zone-cycles" plain.err plain.log
tail -5 dhj.mhd.hst 2>/dev/null
