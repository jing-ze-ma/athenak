#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_1d/reg.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_1d/reg.err
#SBATCH -J rg_reg
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:14:00
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
cd $B/tests_1d/reg_ref
srun -n 1 $B/tests_1d/athena_ref -i $B/inputs/hydro/red_giant_cs.athinput time/nlim=100 problem/ck_data_dir=/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck problem/ck_table=/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck/ck/Premixed_1x_g8_11.txt problem/opac_table=/viper/u2/jinma/ATHENAK/athenak/data/stellar_opac/rosseland_gs98_x0.7_z0.014.txt
cd $B/tests_1d/reg_new
srun -n 1 $B/build_gpu_rg/src/athena -i $B/inputs/hydro/red_giant_cs.athinput time/nlim=100 problem/ck_data_dir=/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck problem/ck_table=/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck/ck/Premixed_1x_g8_11.txt problem/opac_table=/viper/u2/jinma/ATHENAK/athenak/data/stellar_opac/rosseland_gs98_x0.7_z0.014.txt
