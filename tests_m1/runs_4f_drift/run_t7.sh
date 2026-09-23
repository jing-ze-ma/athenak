#!/bin/bash
# usage: run_t7.sh <name> <M5|M7> <N> <be|h2> [overrides]  (one run in drift_0923/<name>)
# The inputs and end states are those of runs_3z_validate (lists/jobs_marshak.txt).
V=/viper/ptmp2/jinma/validate_0923; D=/viper/ptmp2/jinma/drift_0923
X=${X:-/viper/ptmp2/jinma/wt_drift/build_cpu/src/athena}
n=$1; mm=$2; N=$3; sc=$4; shift 4
I=/viper/ptmp2/jinma/wt_drift/tests_m1/runs_4f_drift/inp
inp=$I/radshock_be.athinput; [ $sc = h2 ] && inp=$I/radshock_h2.athinput
if [ $mm = M5 ]; then
  a="time/tlim=2.0e-10 output1/dt=5.0e-11 output2/dt=5.0e-11 problem/m1_shock_rho_r=2.04721116e+01 problem/m1_shock_v_r=2.40600842e+07 problem/m1_shock_t_r=1.86546943e+07 problem/m1_shock_v_l=8.65660331e+07 problem/m1_shock_xs=3.31573616e-02 mesh/x1max=3.35316087e-02 problem/m1_shock_ref=ref_m5.txt rad_m1/implicit_ebath_x1min=1.70873648e+11 rad_m1/implicit_ebath_x1max=9.16223997e+14"
else
  a="time/tlim=4.0e-9 output1/dt=5.0e-10 output2/dt=5.0e-10 problem/m1_shock_rho_r=3.19609192e+01 problem/m1_shock_v_r=2.15758819e+07 problem/m1_shock_t_r=6.75689446e+06 problem/m1_shock_v_l=1.21192446e+08 problem/m1_shock_xs=1.79427684e+00 mesh/x1max=1.82513931e+00 problem/m1_shock_ref=ref_m7p1.txt rad_m1/arad=1.6332399e-42 mesh/x2min=-7.129450e-03 mesh/x2max=7.129450e-03 rad_m1/implicit_ebath_x1min=1.70873648e+15 rad_m1/implicit_ebath_x1max=1.57702166e+17"
fi
d=$D/$n; rm -rf $d; mkdir -p $d; cp $D/ref_*.txt $d/
cd $d && OMP_NUM_THREADS=1 nice -n 10 $X -i $inp $a time/cfl_number=0.4 mesh/nx1=$N \
  meshblock/nx1=$N "$@" > run.log 2>&1; echo "exit=$?" >> run.log
