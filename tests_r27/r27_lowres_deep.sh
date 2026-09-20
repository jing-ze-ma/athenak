#!/bin/bash -l
#SBATCH -J he4deep
#SBATCH -p apu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=06:00:00
#SBATCH --time-min=01:00:00
#
# tests_r27: the DEEP-INNER-WALL low-resolution arm of the 3-D He4 envelope wedge.
#
# Identical to tests_r22/r22_lowres24.sh -- same physics switches, same seed, same
# 96 x 96 angular grid, same chaining -- with ONE change: the inner wall is moved
# from 0.500 R to 0.400 R so the collapsing envelope lands on the deep, thermally
# slow radiative interior instead of on the wall.
#
#   r_in      1.18585e11 (0.500 R)  ->  9.486800e10 (0.400 R)
#   nx1       144                   ->  216   (one MeshBlock in radius; 216+6 <= 520)
#   stretch   refitted by tests_r27/fit2.py on [0.400 R, 1.25 R]: dr = 5.7e8 at the
#             wall, 3.5e8 at 1.02 R, <= 9.2e8 over 0.97-1.10 R, 6.6-13.1e8 in the FeCZ
#             (the 144-cell grid had 6.6-21.1e8 there, so the FeCZ is not coarsened)
#   rad_flux_inner  1.305278e15     ->  2.039496e15  = L_star/(4 pi r_in^2)
#   vdamp_bot_cells 24              ->  74   -- the SAME PHYSICAL TOP (0.7196 R) for the
#             bottom sponge, so the envelope sees exactly the damping it saw before and
#             the new 0.400-0.500 R reservoir is held quiet against the wall
#
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
ARM=$1; shift; A=$B/tests_r27/$ARM; mkdir -p $A; cd $A
OV="problem/inner_bc=wall problem/rt_bottom_flux=true
 problem/mlt_alpha=0.0 problem/rt_top_vacuum=true hydro/vceil_thermalise=true
 problem/vpert_var=eint problem/vpert=3.0e-2 problem/vpert_sp_rand=true
 problem/vdamp_all_mean_only=true problem/vdamp_all_until=18820.0
 problem/vdamp_all_time=500.0
 hydro/eos_rad_rho_hi=1.0e-10 hydro/eos_rad_rho_lo=1.0e-11
 hydro/eos_rad_t_hi=0.0 hydro/eos_rad_t_lo=0.0
 problem/rt_rad_force=false problem/rt_force_tau_gate=false hydro/rad_gate_rho=0.0
 time/cfl_number=0.3"
GRID="problem/ic_profile=${ICF:-$B/tests_r14/ic_he4_tall_neww.txt}
 mesh/x1min=9.486800e10 mesh/x1max=2.964625e11 mesh/nx1=216 meshblock/nx1=216
 mesh/f_stretch_r_c1=-0.391990 mesh/f_stretch_r_c2=+0.936503
 mesh/f_stretch_r_c3=+4.196941 mesh/f_stretch_r_c4=-6.428076
 hydro/rad_flux_inner=2.039496e15 problem/vdamp_bot_cells=74
 mesh/nx2=96 mesh/nx3=96 meshblock/nx2=24 meshblock/nx3=24"
COMMON="output1/dt=47.0 output2/dt=2352.5 output3/dt=940.0 output4/dt=1.0e30
 output5/dt=47.0 problem/rt_profile_dt=47.0 problem/rt_surface_dt=94.0
 problem/column_dump=column_$ARM.txt problem/mlt_dump=mltfaces_$ARM.txt"
if [ -n "$OLDJOB" ] && squeue -h -j $OLDJOB 2>/dev/null | grep -q .; then
  scancel $OLDJOB; sleep 60
fi
RST=$(ls -t $A/rst/*.rst 2>/dev/null | head -1)
NT=$((SLURM_JOB_NUM_NODES*2))
LEFT=$(squeue -h -j $SLURM_JOB_ID -o %L)
SEC=$(echo $LEFT | awk -F'[-:]' '{n=NF; s=$n+60*$(n-1); if(n>2)s+=3600*$(n-2); if(n>3)s+=86400*$(n-3); print s-600}')
if [ "$SEC" -lt 120 ]; then SEC=120; fi
TL=$(printf "%02d:%02d:%02d" $((SEC/3600)) $((SEC%3600/60)) $((SEC%60)))
echo "### nodes $SLURM_JOB_NUM_NODES ranks $NT granted $LEFT athena -t $TL restart ${RST:-none} ic ${ICF:-tall_neww}"
srun -n $NT $B/tests_r11/athena_v24 ${RST:+-r $RST} -i $B/inputs/hydro/he4_presn_sp.athinput \
     -t $TL $COMMON $OV $GRID time/tlim=${TLIM:-47050} "$@"
echo "### athena exit $?"
ls -t $A/rst/*.rst 2>/dev/null | tail -n +2 | xargs -r rm -f
# tlim 47050 s = 10 turnovers.  WATCH: mass fraction inside 0.70 R and inside the
# innermost 10 % of the radius, L_out/L_wall, tot-E drift.
