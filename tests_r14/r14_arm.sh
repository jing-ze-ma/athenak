#!/bin/bash -l
#SBATCH -J r14arm
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
# tests_r14: the t = 0 force check of the radiation redesign.  A copy of
# tests_r13/r13_arm20.sh with rt_profile_dt = 1 s (fbud.py needs the per-second
# profile) and the IC file overridable through ICF.
# usage: [XBIN=...] [ICF=...] sbatch -o <log> r14_arm.sh <armname> [overrides ...]
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
X=${XBIN:-$B/tests_r11/athena_v21}
IC=${ICF:-/viper/u2/jinma/ATHENAK/bench/hestar_presn/ic_he4_presn_sph.txt}
IN=$B/inputs/hydro/he4_presn_sp.athinput
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5"
OUT="output1/dt=10.0 output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30
 output5/dt=1.0e30 problem/rt_profile_dt=1.0 problem/rt_surface_dt=1.0e30
 problem/face_budget=0"
d=$1; shift
A=$B/tests_r14/$d; rm -rf $A; mkdir -p $A; cd $A
echo "=========== ARM $d : X=$X : IC=$IC : $@"
srun -n 1 $X -i $IN -t 00:12:00 $ONED $GNR $OUT time/tlim=300.0 \
  problem/ic_profile=$IC "$@" 2>&1 \
  | grep -E 'cycle=|COLLAPSE|dt is set|FATAL|Terminating|^time=|GATED|TAPERED|floor|c2p|C2P' \
  | tail -60
rm -rf $A/bin $A/rst $A/cbin_hydro_w_2
echo "=========== ARM $d DONE"
