#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r4/g4.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r4/g4.err
#SBATCH -J g4
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
H=$B/tests_r4
COMMON="time/nlim=1 time/tlim=1.0e30 output1/dt=1.0e30 output2/dt=1.0e30
 output3/dt=1.0e30 output4/dt=1.0e30 output5/dt=1.0e30 problem/rt_profile_dt=1.0e30
 problem/rt_surface_dt=1.0e30 problem/mlt_alpha=1.5 problem/vpert=0.0
 problem/mlt_dump=mltfaces.txt problem/column_dump=column.txt"
BASE_D="problem/inner_bc=wall problem/rt_bottom_flux=true  hydro/rad_flux_inner=1.305278e15"
mkdir -p $H/g4_sw && cd $H/g4_sw
srun -n 2 $H/athena_v2 -i $B/inputs/hydro/he4_presn_cs.athinput $COMMON $BASE_D
echo "##### exit $? sweep(new)"
mkdir -p $H/g4_m3 && cd $H/g4_m3
srun -n 2 $H/athena_v2 -i $B/inputs/hydro/he4_presn_cs.athinput $COMMON $BASE_D \
     problem/rt_col3_skip_sweep=true
echo "##### exit $? mode3(new)"
mkdir -p $H/g4_ref && cd $H/g4_ref
srun -n 2 $H/athena_ref -i $B/inputs/hydro/he4_presn_cs.athinput $COMMON $BASE_D
echo "##### exit $? ref"
echo "##### ALL DONE"
