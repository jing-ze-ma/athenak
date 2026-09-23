#!/bin/bash
# run.sh <exe> <outdir> <m0: 2|5> <scheme: be|h2|expl> <N> [extra overrides...]
X=$1; D=$2; M0=$3; S=$4; N=$5; shift 5
I=/viper/ptmp2/jinma/efix_0923/inp
case $S in
  be) IN=$I/radshock_be.athinput ;;
  h2) IN=$I/radshock_h2.athinput ;;
  expl) IN=$I/radshock_be.athinput ;;
esac
M5="problem/m1_shock_rho_r=2.04721116e+01 problem/m1_shock_v_r=2.40600842e+07 problem/m1_shock_t_r=1.86546943e+07 problem/m1_shock_v_l=8.65660331e+07 problem/m1_shock_xs=3.31573616e-02 mesh/x1max=3.35316087e-02 problem/m1_shock_ref=ref_m5.txt"
A="time/tlim=2.0e-10 output1/dt=5.0e-11 output2/dt=5.0e-11 time/cfl_number=0.4 mesh/nx1=$N meshblock/nx1=$N"
[ "$M0" = 5 ] && A="$A $M5"
[ "$S" = expl ] && A="$A rad_m1/transport=explicit"
rm -rf $D; mkdir -p $D; cp $I/ref_m2.txt $I/ref_m5.txt $D/
cd $D && nice -n 10 $X -i $IN -d $D $A "$@" > run.log 2>&1; echo "exit=$?" >> run.log
