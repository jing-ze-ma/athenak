#!/bin/bash -l
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
IN=$R/inputs/tests/rad_m1_advect_pulse.athinput
O=$R/tests_m1/runs_1c
TL="time/tlim=4.8e-6 output1/dt=1.2e-6 output2/dt=1.2e-6"
run () { local n=$1; shift; local d=$O/t4_$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN $TL "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }
# the STATIC counterparts at the SAME opacity, which is what "advected vs static" means
run k500_v0_512   mesh/nx1=512  meshblock/nx1=512  rad_m1/kappa_p=500 &
run k500_v0_1024  mesh/nx1=1024 meshblock/nx1=1024 rad_m1/kappa_p=500 &
# the "static" QUOKKA case is kappa = 100 with beta = 3.3e-5 (v = 1e6)
run k100_v1e6_512 mesh/nx1=512  meshblock/nx1=512  problem/pulse_v=1.0e6 &
run k100_v1e6_1024 mesh/nx1=1024 meshblock/nx1=1024 problem/pulse_v=1.0e6 &
run k100_v0_1024b mesh/nx1=1024 meshblock/nx1=1024 &
run d512_none_nosplit mesh/nx1=512 meshblock/nx1=512 rad_m1/kappa_p=500 \
    problem/pulse_v=3.0e7 rad_m1/thick_flux=none rad_m1/advect_split=false &
wait
echo T4X_DONE
