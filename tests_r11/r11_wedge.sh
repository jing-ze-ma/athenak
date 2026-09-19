#!/bin/bash -l
#SBATCH -J he4wedge
#SBATCH -p apu
#SBATCH --nodes=1-4
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=06:00:00
#SBATCH --time-min=01:00:00
# tests_r11: the 90x90 deg spherical-polar wedge, 96x256x256, recipe of w3d, toward 5
# turnovers.  FLEXIBLE for backfill on a full partition: 1-4 nodes, and the scheduler may
# shorten the limit down to 1 h; athena's own -t is taken from what was actually granted.
# Resubmitting the same command chains from the newest restart file.
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
# UPGRADE: sbatch --export=ALL,OLDJOB=<id> --nodes=2-4 ... takes over from a smaller running
# job of the same arm: cancel it, let its files settle, continue from the newest restart.
if [ -n "$OLDJOB" ] && squeue -h -j $OLDJOB 2>/dev/null | grep -q .; then
  scancel $OLDJOB; sleep 60
fi
RST=$(ls -t $A/rst/*.rst 2>/dev/null | head -1)
NT=$((SLURM_JOB_NUM_NODES*2))
LEFT=$(squeue -h -j $SLURM_JOB_ID -o %L)              # [d-]hh:mm:ss or mm:ss granted
SEC=$(echo $LEFT | awk -F'[-:]' '{n=NF; s=$n+60*$(n-1); if(n>2)s+=3600*$(n-2); if(n>3)s+=86400*$(n-3); print s-600}')
TL=$(printf "%02d:%02d:%02d" $((SEC/3600)) $((SEC%3600/60)) $((SEC%60)))
echo "### nodes $SLURM_JOB_NUM_NODES ranks $NT granted $LEFT athena -t $TL restart ${RST:-none}"
srun -n $NT $B/tests_r11/athena_v18 ${RST:+-r $RST} -i $B/inputs/hydro/he4_presn_sp.athinput \
     -t $TL $COMMON $OV time/tlim=23525 "$@"
echo "### athena exit $?"
ls -t $A/rst/*.rst 2>/dev/null | tail -n +2 | xargs -r rm -f
# the 3-D dumps (bin every 0.5 turnover) are KEPT: they are the convection diagnostics
