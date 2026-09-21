#!/bin/bash -l
# The runs behind tests_m1/plots/*.json.  Default scheme unless stated:
# thick_flux = ap_hll, ap_form = alpha2, reconstruct = plm, advect_split = true.
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
O=$R/tests_m1/runs_1cB
IT=$R/inputs/tests/rad_m1_thick_pulse.athinput
IJ=$R/inputs/tests/rad_m1_jump.athinput
IM=$R/inputs/tests/rad_m1_marshak.athinput
run () { local n=$1; local IN=$2; shift 2; local d=$O/$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }

# (iii) T3 thick pulse at tau_cell = 1e3, the three thick_flux choices
for tf in ap_hll scaled none; do
  run pl_t3_$tf $IT rad_m1/kappa_s=128000 rad_m1/thick_flux=$tf time/tlim=12.0 \
      output1/dt=12.0 output2/dt=12.0 &
done
# (iv) T3c clipped top-hat at tau_cell = 1e4
for tf in ap_hll none; do
  run pl_tophat_$tf $IT rad_m1/kappa_s=1280000 rad_m1/thick_flux=$tf \
      problem/m1_test=tophat time/tlim=100.0 output1/dt=100.0 output2/dt=100.0 &
done
# (vi) T6 Marshak at 64 cells
for tf in ap_hll none; do
  run pl_marshak_$tf $IM mesh/nx1=64 meshblock/nx1=64 rad_m1/thick_flux=$tf &
done
# (vii) T3b case A, tau_cell 1 -> 1e3
run pl_t3b_ap_hll $IJ rad_m1/kappa_f=64.0 problem/jump_flux=1.0e-6 time/tlim=1.0e5 \
    output1/dt=5.0e4 output2/dt=5.0e4 &
wait
echo PLOTS_RUNS_DONE
