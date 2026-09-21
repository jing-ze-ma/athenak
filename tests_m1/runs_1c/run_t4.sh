#!/bin/bash -l
# T4 (design sect. 8): the advected thick pulse, static and dynamic.
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
IN=$R/inputs/tests/rad_m1_advect_pulse.athinput
O=$R/tests_m1/runs_1c
TL=4.8e-6
DT=1.2e-6

run () {
  local n=$1; shift
  local d=$O/t4_$n
  rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN time/tlim=$TL output1/dt=$DT output2/dt=$DT "$@" \
      > run.log 2>&1; echo "exit=$?" >> run.log)
}

DYN="rad_m1/kappa_p=500 problem/pulse_v=3.0e7"
N512="mesh/nx1=512 meshblock/nx1=512"
N1024="mesh/nx1=1024 meshblock/nx1=1024"

case $1 in
a)
  run s512_aphll_split        $N512
  run d512_aphll_split        $N512 $DYN
  run d512_aphll_nosplit      $N512 $DYN rad_m1/advect_split=false
  ;;
b)
  run d512_scaled_split       $N512 $DYN rad_m1/thick_flux=scaled
  run d512_none_split         $N512 $DYN rad_m1/thick_flux=none
  run d512_aphll_split_ovc    $N512 $DYN rad_m1/source_form=ovc
  run d512_aphll_split_nsub1  $N512 $DYN rad_m1/subcycle=false
  ;;
c)
  run s1024_aphll_split       $N1024
  run d1024_aphll_split       $N1024 $DYN
  run d1024_aphll_nosplit     $N1024 $DYN rad_m1/advect_split=false
  ;;
esac
echo T4_DONE_$1
