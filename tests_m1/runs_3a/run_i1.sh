#!/bin/bash -l
# I1: the thick pulse with IMPLICIT transport, at tau_cell = 10, 1e3, 1e6 and a scan of
# the radiation CFL c dt/dx, plus the Nyquist mode.  Each case is run, analysed and its
# dumps deleted immediately (the viper inode quota is at its limit).
# <time>/cfl_number is raised so the HYDRO step never limits the radiation step; the gas
# is uniform and static (gas_feedback = false), so that is inert.
R=/viper/u2/jinma/ATHENAK/bench/wt_m1impl
X=$R/build_cpu_m1/src/athena
T=$R/tests_m1
O=$T/tests_out_i1
IT=$R/inputs/tests/rad_m1_thick_pulse.athinput
mkdir -p $O
cd $T

one () {  # name kappa tlim odt extra...
  local n=$1 kap=$2 tl=$3 odt=$4; shift 4
  rm -rf $O/$n; mkdir -p $O/$n
  (cd $O/$n && $X -i $IT time/cfl_number=1e6 rad_m1/transport=implicit_x1 \
      rad_m1/kappa_s=$kap time/tlim=$tl output1/dt=$odt output2/dt=1.0e30 "$@" \
      > run.log 2>&1)
  local ex=$?
  local it=$(grep -o "iterations mean=[0-9e.+-]*" $O/$n/run.log | tail -1)
  local nc=$(grep -o "NON-CONVERGED=[0-9e.+-]*" $O/$n/run.log | tail -1)
  local nstep=$(grep -o "solves=[0-9e.+-]*" $O/$n/run.log | tail -1)
  printf "%-22s exit=%s %s %s %s  " $n $ex "$nstep" "$it" "$nc"
  shift 0
  python3 $NYQ t3_pulse.py $O/$n/bin/*.bin --c 1 --rho 1 --kappa $kap \
      --subtract-min --quiet $NYQFLAG
  rm -rf $O/$n/bin $O/$n/tab
}

for cfl in 0.4 1 10 100 10000; do
  NYQFLAG=""
  one i1_tau1e1_c$cfl 1280      0.5    0.05  rad_m1/implicit_cfl=$cfl
  one i1_tau1e3_c$cfl 128000    12.0   1.2   rad_m1/implicit_cfl=$cfl
  one i1_tau1e6_c$cfl 128000000 3000.0 300.0 rad_m1/implicit_cfl=$cfl
  NYQFLAG="--nyquist"
  one i1_nyq_c$cfl    128000    20.0   2.0   rad_m1/implicit_cfl=$cfl \
      problem/nyquist_amp=1.0e-3
done
echo I1_DONE
