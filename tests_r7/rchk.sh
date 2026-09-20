#!/bin/bash -l
#SBATCH -J r7_rchk
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:14:00
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/rchk.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/rchk.err
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5"
OUT="output1/dt=23.5 output2/dt=1.0e30 output4/dt=1.0e30 output5/dt=1.0e30
 problem/rt_profile_dt=1.0e30 problem/rt_surface_dt=1.0e30 problem/face_budget=0"
# The chain, with the OLD binary (no closure state in the file) and with the NEW one.
# rt_rad_force=false: tests_1d D(iv) showed the chained restart is reproducible to 1 ulp
# in that configuration, so any larger difference here is the closure re-seed and nothing
# else.  NF selects it.
chain () {
  tag=$1; X=$2; shift 2
  A=$B/tests_r7/rc_${tag}_cont; rm -rf $A; mkdir -p $A; cd $A
  echo "########## $tag CONT ($X) $@"
  srun -n 1 $X -i $B/inputs/hydro/he4_presn_cs.athinput -t 00:02:30 $ONED $GNR $OUT \
    output3/dt=2352.5 time/tlim=4705.0 problem/column_dump=col.txt \
    problem/mlt_dump=mlt_${tag}_cont.txt "$@" 2>&1 | grep -E 'MLT closure|Terminating'
  A2=$B/tests_r7/rc_${tag}_rst; rm -rf $A2; mkdir -p $A2; cd $A2
  echo "########## $tag RST ($X)"
  srun -n 1 $X -r $A/rst/he4.00001.rst -i $B/inputs/hydro/he4_presn_cs.athinput \
    -t 00:02:30 $ONED $GNR $OUT output3/dt=1.0e30 time/tlim=4705.0 \
    problem/column_dump=col.txt problem/mlt_dump=mlt_${tag}_rst.txt "$@" 2>&1 \
    | grep -E 'MLT closure|Terminating'
  echo "---------- $tag hst"
  python3 $B/tests_r7/cmphst.py $A/he4.hydro.hst $A2/he4.hydro.hst
  python3 $B/tests_r7/firstline.py $A/he4.hydro.hst $A2/he4.hydro.hst
  rm -rf $A/bin $A2/bin $A/rst $A2/rst $A/cbin_hydro_w_2 $A2/cbin_hydro_w_2
}
chain oldNF $B/tests_r4/athena_v2 problem/rt_rad_force=false
chain newNF $B/tests_r7/athena_v3 problem/rt_rad_force=false
chain oldF  $B/tests_r4/athena_v2
chain newF  $B/tests_r7/athena_v3
