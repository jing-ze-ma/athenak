#!/bin/bash -l
#SBATCH -J r7_chk
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:14:00
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/chk.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/chk.err
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
X=$B/tests_r7/athena_v3
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5"
OUT="output1/dt=23.5 output2/dt=1.0e30 output4/dt=1.0e30 output5/dt=1.0e30
 problem/rt_profile_dt=1.0e30 problem/rt_surface_dt=1.0e30 problem/face_budget=0"
# ---------------- B.  the restart chain: 1.0 turnover in one piece vs 0.5 + 0.5
A=$B/tests_r7/rb_cont; rm -rf $A; mkdir -p $A; cd $A
echo "=========== rb_cont: one piece to 1.0 turnover, rst written at 0.5"
srun -n 1 $X -i $B/inputs/hydro/he4_presn_cs.athinput -t 00:04:00 $ONED $GNR $OUT \
  output3/dt=2352.5 time/tlim=4705.0 problem/column_dump=col.txt \
  problem/mlt_dump=mlt_cont.txt 2>&1 | tail -4
ls -la $A/rst/
A2=$B/tests_r7/rb_rst; rm -rf $A2; mkdir -p $A2; cd $A2
R=$(ls $A/rst/*00001.rst)
echo "=========== rb_rst: restart from $R to 1.0 turnover"
srun -n 1 $X -r $R -i $B/inputs/hydro/he4_presn_cs.athinput -t 00:04:00 \
  $ONED $GNR $OUT output3/dt=1.0e30 time/tlim=4705.0 \
  problem/column_dump=col.txt problem/mlt_dump=mlt_rst.txt 2>&1 | tail -6
echo "=========== the hst comparison"
python3 $B/tests_r7/cmphst.py $A/he4.hydro.hst $A2/he4.hydro.hst
rm -rf $A/bin $A2/bin $A/cbin_hydro_w_2 $A2/cbin_hydro_w_2
# ---------------- A.  the two survivors, out to 2.5 turnovers
run () {
  d=$1; shift
  AA=$B/tests_r7/$d; mkdir -p $AA; cd $AA
  echo "=========== ARM $d : $@"
  srun -n 1 $X -i $B/inputs/hydro/he4_presn_cs.athinput -t 00:03:00 $ONED $GNR $OUT \
    output3/dt=1.0e30 time/tlim=11762.0 problem/column_dump=column_${d}.txt \
    problem/mlt_dump=mltfaces_${d}.txt "$@" 2>&1 \
    | grep -E 'cycle=|COLLAPSE|dt is set|^    r=|FATAL|Terminating|^time=' | tail -20
  rm -rf $AA/bin $AA/rst $AA/cbin_hydro_w_2
}
run x_noforce   problem/rt_rad_force=false
run x_nostrang  problem/rt_strang=false
run x_base
