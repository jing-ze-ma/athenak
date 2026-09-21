#!/bin/bash -l
# T4b (design sect. 8): uniform advected medium, beta tau_cell = 1e3.
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
IN=$R/inputs/tests/rad_m1_advect_uniform.athinput
O=$R/tests_m1/runs_1c

run () {
  local n=$1; shift
  local d=$O/t4b_$n
  rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log)
}

run aphll_split_full  rad_m1/advect_split=true  rad_m1/source_form=full
run aphll_split_ovc   rad_m1/advect_split=true  rad_m1/source_form=ovc
run aphll_nosplit_full rad_m1/advect_split=false rad_m1/source_form=full
run aphll_nosplit_ovc rad_m1/advect_split=false rad_m1/source_form=ovc
run none_split_full   rad_m1/thick_flux=none    rad_m1/advect_split=true
run scaled_split_full rad_m1/thick_flux=scaled  rad_m1/advect_split=true
echo T4B_DONE
