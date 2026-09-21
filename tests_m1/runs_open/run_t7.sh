#!/bin/bash -l
# T7 (design sect. 8): Lowrie & Edwards radiative shocks, Mach 2 subcritical and Mach 5
# supercritical, started FROM the semi-analytic profile and held in the shock frame by
# Dirichlet ends.  Gate: L1 < 2 % in rho, T_gas, T_rad.  Controls: closure = eddington
# (Skinner & Ostriker found the diffusion closure more accurate in these thick shocks)
# and chat/c = 0.1, 0.01 (the RSLA error).
R=/viper/u2/jinma/ATHENAK/bench/wt_m1gates
X=$R/build/src/athena
IN=$R/inputs/tests/rad_m1_radshock.athinput
T=$R/tests_m1
O=$T/runs_open
M5="problem/m1_shock_rho_r=2.04721116e+01 problem/m1_shock_v_r=2.40600842e+07 \
problem/m1_shock_t_r=1.86546943e+07 problem/m1_shock_v_l=8.65660331e+07 \
problem/m1_shock_xs=3.31573616e-02 mesh/x1max=3.35316087e-02 \
problem/m1_shock_ref=ref_m5.txt"
TL="time/tlim=1.0e-10 output1/dt=5.0e-11 output2/dt=5.0e-11"
(cd $T && python3 t7_radshock.py --m0 2 --write-ref $O/ref_m2.txt --quiet > /dev/null)
(cd $T && python3 t7_radshock.py --m0 5 --write-ref $O/ref_m5.txt --quiet > /dev/null)
run () { local n=$1; shift; local d=$O/t7_$n; rm -rf $d; mkdir -p $d
  cp $O/ref_m2.txt $O/ref_m5.txt $d/
  (cd $d && $X -i $IN $TL "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
N5="mesh/nx1=1024 meshblock/nx1=1024"
run m2_c1       &
run m2_edd      rad_m1/closure=eddington &
run m2_chat0.1  rad_m1/chat_over_c=0.1 &
run m2_chat0.01 rad_m1/chat_over_c=0.01 &
run m2_sub      rad_m1/subcycle=true &
run m5_c1       $M5 mesh/nx1=2048 meshblock/nx1=2048 &
run m5_c1_n1024 $M5 $N5 &
run m5_edd      $M5 $N5 rad_m1/closure=eddington &
run m5_chat0.1  $M5 $N5 rad_m1/chat_over_c=0.1 &
run m5_chat0.01 $M5 $N5 rad_m1/chat_over_c=0.01 &
wait
echo T7_DONE
