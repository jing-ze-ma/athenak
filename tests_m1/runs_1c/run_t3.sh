#!/bin/bash -l
# T3 family (design sect. 8 / 10): the AP gate, for both ap_form choices.
# Usage: run_t3.sh <ap_form>            (unified | alpha2)
A=$1
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
IN=$R/inputs/tests/rad_m1_thick_pulse.athinput
O=$R/tests_m1/runs_1c

run () {  # name  extra...
  local n=$1; shift
  local d=$O/t3_${A}_${n}
  rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN rad_m1/ap_form=$A "$@" > run.log 2>&1; echo "exit=$?" >> run.log)
}

for rec in dc plm; do
  run tau1e1_$rec rad_m1/kappa_s=1280     rad_m1/reconstruct=$rec time/tlim=0.5 \
      output1/dt=0.05
  run tau1e3_$rec rad_m1/kappa_s=128000   rad_m1/reconstruct=$rec time/tlim=12.0 \
      output1/dt=1.2
  run tau1e6_$rec rad_m1/kappa_s=128000000 rad_m1/reconstruct=$rec time/tlim=3000.0 \
      output1/dt=300.0
done
# resolution order at tau_cell = 10 (plm): 256 cells at the same tau_cell
run tau1e1_plm_n256 rad_m1/kappa_s=2560 mesh/nx1=256 meshblock/nx1=256 \
    time/tlim=0.5 output1/dt=0.05
# Nyquist mode at tau_cell = 1e3
run nyq_plm rad_m1/kappa_s=128000 problem/nyquist_amp=1.0e-3 time/tlim=20.0 \
    output1/dt=2.0
# T3c clipped top-hat at tau_cell = 1e4
run tophat rad_m1/kappa_s=1280000 problem/m1_test=tophat time/tlim=100.0 \
    output1/dt=10.0
echo "T3_DONE_$A"
