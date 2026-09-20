#!/bin/bash -l
#SBATCH -J he4lowres
#SBATCH -p apu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=06:00:00
#SBATCH --time-min=01:00:00
#
# tests_r22: the LOW-RESOLUTION RELAXATION arm of the 3-D He4 envelope wedge.
#
# Same physics, same IC, same seed as the production wedge (tests_r11/r11_wedge23.sh +
# the argument list of job 11865437 = wedge11) -- ONLY the angular resolution changes:
#
#     production   144 x 192 x 192   blocks 144 x 24 x 24   (8x8 = 64 blocks)   5.31e6 cells
#     this arm     144 x  96 x  96   blocks 144 x 24 x 24   (4x4 = 16 blocks)   1.33e6 cells
#
# WHY 96 AND NOT 64 (see tests_r22/README.md for the full argument):
#   * 90 deg / 192 = 0.469 deg = 1.94e9 cm per cell at R = 2.3717e11 cm.
#     90 deg /  96 = 0.938 deg = 3.88e9 cm;   90 deg / 64 = 1.41 deg = 5.82e9 cm.
#   * the ENTROPY SEED (vpert_sp_rand, vpert_sp_kmax = 8) has a shortest wavelength of
#     90/8 = 11.25 deg = 12 cells at 96 and 8 cells at 64: both grids carry the seed.
#   * the POROUS CHANNELS measured at 3.5 turnovers have autocorrelation lengths
#     5-14e9 cm = 0.5-1.2 deg: 3-9 cells at 192, 1.3-3.6 cells at 96, 0.9-2.4 at 64.
#     96 therefore still carries the LARGEST channels (marginally, 3-4 cells across) and
#     nothing of the small ones; 64 carries no channel at all -- it would relax the mean
#     state under a convective flux that the grid cannot produce.  96 is the coarsest
#     grid on which the relaxation is still driven by (under-resolved) convection rather
#     than by numerical diffusion alone.
#   * 4x4 angular blocks is well inside the implicit angular ADI limit NADIB = 8, and
#     144 radial cells stay in a single MeshBlock as the column solve requires.
#     16 blocks / 2 GPUs = 8 blocks per GPU (production: 64/8 = 8 -- identical per-GPU
#     block count, so the per-GPU work is 1/1 in blocks and the arm is ~4x cheaper in
#     total cells while using 1/4 of the nodes).
#
# ONE NODE, 2 GPUs, chained exactly like r11_wedge23.sh: resubmitting the same command
# line continues from the arm's own newest restart file.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
ARM=$1; shift; A=$B/tests_r22/$ARM; mkdir -p $A; cd $A
# the production physics switches (r11_wedge23.sh's OV, with the wedge11 command line's
# mlt_alpha = 0 and the entropy seed folded in)
OV="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=0.0 problem/rt_top_vacuum=true
 problem/vpert_var=eint problem/vpert=3.0e-2 problem/vpert_sp_rand=true
 problem/vdamp_all_mean_only=true problem/vdamp_all_until=18820.0
 problem/vdamp_all_time=500.0
 hydro/eos_rad_rho_hi=1.0e-10 hydro/eos_rad_rho_lo=1.0e-11
 hydro/eos_rad_t_hi=0.0 hydro/eos_rad_t_lo=0.0
 problem/rt_rad_force=false problem/rt_force_tau_gate=false hydro/rad_gate_rho=0.0
 time/cfl_number=0.3"
# THE GRID: the tall IC, 144 radial cells in one block, 96x96 angular in 24x24 blocks
GRID="problem/ic_profile=$B/tests_r14/ic_he4_tall_neww.txt
 mesh/x1max=2.964625e11 mesh/nx1=144 meshblock/nx1=144
 mesh/f_stretch_r_c1=0.992525 mesh/f_stretch_r_c2=-0.404821
 mesh/f_stretch_r_c3=-0.372255 mesh/f_stretch_r_c4=-1.499759
 mesh/nx2=96 mesh/nx3=96 meshblock/nx2=24 meshblock/nx3=24"
COMMON="output1/dt=47.0 output2/dt=2352.5 output3/dt=940.0 output4/dt=1.0e30
 output5/dt=47.0 problem/rt_profile_dt=47.0 problem/rt_surface_dt=94.0
 problem/column_dump=column_$ARM.txt problem/mlt_dump=mltfaces_$ARM.txt"
if [ -n "$OLDJOB" ] && squeue -h -j $OLDJOB 2>/dev/null | grep -q .; then
  scancel $OLDJOB; sleep 60
fi
RST=$(ls -t $A/rst/*.rst 2>/dev/null | head -1)
NT=$((SLURM_JOB_NUM_NODES*2))
LEFT=$(squeue -h -j $SLURM_JOB_ID -o %L)              # [d-]hh:mm:ss or mm:ss granted
SEC=$(echo $LEFT | awk -F'[-:]' '{n=NF; s=$n+60*$(n-1); if(n>2)s+=3600*$(n-2); if(n>3)s+=86400*$(n-3); print s-600}')
if [ "$SEC" -lt 120 ]; then SEC=120; fi
TL=$(printf "%02d:%02d:%02d" $((SEC/3600)) $((SEC%3600/60)) $((SEC%60)))
echo "### nodes $SLURM_JOB_NUM_NODES ranks $NT granted $LEFT athena -t $TL restart ${RST:-none}"
srun -n $NT $B/tests_r11/athena_v23 ${RST:+-r $RST} -i $B/inputs/hydro/he4_presn_sp.athinput \
     -t $TL $COMMON $OV $GRID time/tlim=47050 "$@"
echo "### athena exit $?"
ls -t $A/rst/*.rst 2>/dev/null | tail -n +2 | xargs -r rm -f
# tlim 47050 s = 10 turnovers: the relaxation target.  Stop earlier by the criteria in
# tests_r22/README.md (L_out/L within a few % of 1, tot-E flat over >= 1 turnover).
