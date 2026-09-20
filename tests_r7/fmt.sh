#!/bin/bash -l
#SBATCH -J r7_fmt
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:10:00
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/fmt.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/fmt.err
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
C="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 output1/dt=1.0e30 output2/dt=1.0e30 output4/dt=1.0e30 output5/dt=1.0e30
 problem/rt_profile_dt=1.0e30 problem/rt_surface_dt=1.0e30 problem/face_budget=0
 output3/dt=1.0 time/nlim=5 time/tlim=1.0e30"
# THE FORMAT GATE: with problem/mlt_alpha = 0 the new binary must write the SAME BYTES the
# old one does, and each binary must read the other's file.
for tag in old new; do
  [ $tag = old ] && X=$B/tests_r4/athena_v2 || X=$B/tests_r7/athena_v3
  A=$B/tests_r7/fmt_$tag; rm -rf $A; mkdir -p $A; cd $A
  srun -n 1 $X -i $B/inputs/hydro/he4_presn_cs.athinput $ONED $C problem/mlt_alpha=0.0 \
    problem/column_dump=c.txt > run.log 2>&1
  echo "$tag mlt_alpha=0 rst: $(ls -la $A/rst/)"
done
echo "=== cmp of the two mlt_alpha=0 restart files (must be identical) ==="
cmp $B/tests_r7/fmt_old/rst/he4.00001.rst $B/tests_r7/fmt_new/rst/he4.00001.rst \
  && echo "BITWISE IDENTICAL"
echo "=== the new binary reads the OLD binary's file (mlt_alpha=0) ==="
cd $B/tests_r7/fmt_new && srun -n 1 $B/tests_r7/athena_v3 \
  -r $B/tests_r7/fmt_old/rst/he4.00001.rst -i $B/inputs/hydro/he4_presn_cs.athinput \
  $ONED $C problem/mlt_alpha=0.0 time/nlim=8 2>&1 | grep -E "MLT closure|cycle=|Termin"
echo "=== the new binary reads the OLD binary's file WITH mlt_alpha=1.5 (must warn) ==="
srun -n 1 $B/tests_r7/athena_v3 -r $B/tests_r7/fmt_old/rst/he4.00001.rst \
  -i $B/inputs/hydro/he4_presn_cs.athinput $ONED $C problem/mlt_alpha=1.5 time/nlim=8 \
  2>&1 | grep -E "MLT closure|cycle=|Termin"
echo "=== the OLD binary reads a NEW mlt_alpha=1.5 file (must FAIL cleanly or be refused) ==="
A=$B/tests_r7/fmt_newm; rm -rf $A; mkdir -p $A; cd $A
srun -n 1 $B/tests_r7/athena_v3 -i $B/inputs/hydro/he4_presn_cs.athinput $ONED $C \
  problem/mlt_alpha=1.5 problem/column_dump=c.txt problem/mlt_dump=m.txt > run.log 2>&1
srun -n 1 $B/tests_r4/athena_v2 -r $A/rst/he4.00001.rst \
  -i $B/inputs/hydro/he4_presn_cs.athinput $ONED $C problem/mlt_alpha=1.5 time/nlim=8 \
  2>&1 | grep -E "FATAL|broken|cycle=|Termin" | head -5
rm -rf $B/tests_r7/fmt_old/bin $B/tests_r7/fmt_new/bin $A/bin
