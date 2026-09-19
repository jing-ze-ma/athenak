#!/bin/bash -l
#SBATCH -J r11_3da
#SBATCH -p apu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=03:30:00
# tests_r11: the 3-D smoke arm (tests_r8 Gnr8 configuration: mlt_alpha 1.5, wall + bottom
# flux, T-gated taper, floors, vpert 1e-3) run in ONE piece toward 5 turnovers with the
# current binary (athena_v15 = a6d66c3c + rt_top_vacuum plumbing), restart every half
# turnover so it can be chained.  Question: does resolved convection start and does the
# envelope settle, once nothing kills the run?
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
ARM=$1; shift; A=$B/tests_r11/$ARM; mkdir -p $A; cd $A
OV="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5 problem/vpert=1.0e-3"
COMMON="output1/dt=47.0 output2/dt=2352.5 output3/dt=940.0 output4/dt=1.0e30
 output5/dt=47.0 problem/rt_profile_dt=47.0 problem/rt_surface_dt=94.0
 problem/column_dump=column_$ARM.txt problem/mlt_dump=mltfaces_$ARM.txt"
RST=$(ls -t $A/rst/*.rst 2>/dev/null | head -1)
if [ -n "$RST" ]; then
  srun -n 2 $B/tests_r11/athena_v15 -r $RST -i $B/inputs/hydro/he4_presn_cs.athinput \
       -t 03:20:00 $COMMON $OV time/tlim=9410.0 "$@"
else
  srun -n 2 $B/tests_r11/athena_v15 -i $B/inputs/hydro/he4_presn_cs.athinput \
       -t 03:20:00 $COMMON $OV time/tlim=9410.0 "$@"
fi
echo "### athena exit $?"
ls -t $A/rst/*.rst 2>/dev/null | tail -n +2 | xargs -r rm -f
rm -rf $A/bin $A/cbin
