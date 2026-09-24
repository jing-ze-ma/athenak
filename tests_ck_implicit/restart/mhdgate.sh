#!/bin/bash
# MHD restart gate (CPU, 12 ranks, login node nice 10): the prod4 restart (read in place)
# with c2 + ck_impl_every = 4: 6 cycles straight vs 3 + restart + 3.
# usage: mhdgate.sh <bin> <rst>
Q=/viper/ptmp2/jinma/ckrst_0924
T=/viper/ptmp2/jinma/wt_ckrst/tests_ck_implicit/restart
B=$1; RST=$2
T4="problem/ck_implicit=true problem/ck_impl_arat=1e30 problem/ck_impl_frozen_op=true problem/ck_impl_lin=true problem/ck_impl_lin_thr=1 problem/ck_impl_fuse=true problem/ck_impl_jac_lin=true problem/ck_impl_cvsec=true problem/ck_impl_verbose=true"
C2="$T4 problem/ck_impl_once=true problem/ck_impl_xstep=2 problem/ck_impl_jreuse=0.2 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 problem/ck_impl_nosync=true"
X="$C2 problem/ck_impl_every=4"
O="output1/dt=1e-6 output2/dt=1e30 output3/dt=1e-6 output4/dt=1e30 time/tlim=1e30"
module purge; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
export OMP_NUM_THREADS=1
read T0 DT0 NC0 <<< $(python3 $T/rstinfo.py $RST)
G=$Q/mhdgate
for r in s6 a3; do d=$G/$r; rm -rf $d; mkdir -p $d; cd $d
  nice -n 10 mpirun -np 12 $B -r $RST -i $T/prod4_mhd_rst.athinput $O $X time/nlim=$((NC0 + ${r:1})) > run.log 2>&1 &
done
wait
d=$G/r3; rm -rf $d; mkdir -p $d; cd $d
nice -n 10 mpirun -np 12 $B -r $(ls $G/a3/rst/*.rst | tail -1) $O time/nlim=$((NC0 + 6)) > run.log 2>&1
echo "RST $RST t=$T0 ncycle=$NC0" > $G/gate.out
python3 $T/cmp.py rst $G/s6 $G/r3 >> $G/gate.out
grep -h "ck restart\|implicit-ck" $G/r3/run.log >> $G/gate.out
echo MHDGATE_DONE >> $G/gate.out
