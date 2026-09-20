#!/bin/bash -l
#SBATCH -J he4_budget
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
H=$B/tests_3d/handover

COMMON="time/nlim=1 time/tlim=1.0e30 output1/dt=1.0e30 output2/dt=1.0e30
 output3/dt=1.0e30 output4/dt=1.0e30 output5/dt=1.0e30 problem/rt_profile_dt=1.0e30
 problem/rt_surface_dt=1.0e30 problem/mlt_alpha=1.5 problem/vpert=0.0
 problem/mlt_dump=mltfaces.txt problem/column_dump=column.txt"

BASE_H="problem/inner_bc=wall problem/rt_bottom_flux=false hydro/rad_flux_inner=1.305278e15"
BASE_D="problem/inner_bc=wall problem/rt_bottom_flux=true  hydro/rad_flux_inner=1.305278e15"

run () {
  local NAME="$1"; shift
  local OV="$*"
  mkdir -p $H/$NAME && cd $H/$NAME
  echo "##### RUN $NAME  OV: $OV"
  srun -n 2 $B/build_gpu_rg/src/athena -i $B/inputs/hydro/he4_presn_cs.athinput \
       $COMMON $OV
  echo "##### exit $? for $NAME"
}

run H1 "$BASE_H hydro/rad_tau_lo=20   hydro/rad_tau_hi=300"
run H2 "$BASE_H hydro/rad_tau_lo=5    hydro/rad_tau_hi=50"
run H3 "$BASE_H hydro/rad_tau_lo=50   hydro/rad_tau_hi=1000"
run M0 "$BASE_D problem/rt_impl_mixed=0"
run Q4 "$BASE_D problem/ck_nquad=4"
echo "##### ALL DONE"
