#!/bin/bash -l
#SBATCH -J r7_fmt2
#SBATCH -p apudev
#SBATCH --nodes=1 --ntasks=1 --gres=gpu:1 --cpus-per-task=8 --time=00:08:00
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/fmt2.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/fmt2.err
module purge; module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
C="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 output1/dt=1.0e30 output2/dt=1.0e30 output4/dt=1.0e30 output5/dt=1.0e30
 problem/rt_profile_dt=1.0e30 problem/rt_surface_dt=1.0e30 problem/face_budget=0
 output3/dt=1.0 time/nlim=5 time/tlim=1.0e30 problem/mlt_alpha=0.0"
# CONTROL: the SAME binary twice.  If its own two files differ in the same 16 bytes, the
# difference the old-vs-new comparison found is uninitialised struct padding in the header
# and not this change.
for r in a b; do
  A=$B/tests_r7/fmt_ctl$r; rm -rf $A; mkdir -p $A; cd $A
  srun -n 1 $B/tests_r4/athena_v2 -i $B/inputs/hydro/he4_presn_cs.athinput $ONED $C \
    problem/column_dump=c.txt > run.log 2>&1
done
echo "=== athena_v2 against ITSELF, two runs ==="
cmp $B/tests_r7/fmt_ctla/rst/he4.00001.rst $B/tests_r7/fmt_ctlb/rst/he4.00001.rst \
  && echo "BITWISE IDENTICAL"
python3 - <<'PY'
import sys
a=open('/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/fmt_ctla/rst/he4.00001.rst','rb').read()
b=open('/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/fmt_ctlb/rst/he4.00001.rst','rb').read()
d=[i for i in range(min(len(a),len(b))) if a[i]!=b[i]]
print('sizes',len(a),len(b),'ndiff',len(d),'range',(d[0],d[-1]) if d else None)
PY
rm -rf $B/tests_r7/fmt_ctla/bin $B/tests_r7/fmt_ctlb/bin
