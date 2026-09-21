#!/bin/bash -l
# Milestone 1d: the 1c gate table for the DEFAULT f_source = cell and for f_source = wb.
cd "$(dirname "$0")/.."
R=runs_open/reg
A5="--c 2.99792458e10 --rho 1.2 --kappa 500 --v 3.0e7 --drift-tol 1 --quiet"
T6="--table runs_open/t6_ref_sn.txt --marshak --a-rad 1e30 --cv 1.5 --rho 1 --centre 0.0 --quiet"
for F in cell wb; do
  for n in tau1e1:1280 tau1e3:128000 tau1e6:128000000; do
    k=${n#*:}; printf "%-5s t3_%-9s " $F ${n%:*}
    python3 t3_pulse.py $R/${F}_t3_${n%:*}/bin/*.bin --c 1 --rho 1 --kappa $k \
        --subtract-min --quiet
  done
  printf "%-5s t3_nyquist   " $F
  python3 t3_pulse.py $R/${F}_t3_nyq/bin/*.bin --c 1 --rho 1 --kappa 128000 \
      --nyquist --quiet
  printf "%-5s t3c_tophat   " $F
  python3 t3_pulse.py $R/${F}_tophat/bin/*.bin --c 1 --rho 1 --kappa 1280000 \
      --subtract-min --quiet
  S=$R/${F}_t4_s512/tab/m1_advect_pulse.m1.00004.tab
  for n in t4_d512 t4_d512ns; do
    d=$R/${F}_$n; printf "%-5s %-12s " $F $n
    python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab $A5 --static $S
  done
  for n in t4b_full t4b_ovc; do
    d=$R/${F}_$n; printf "%-5s %-12s " $F $n
    python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab \
        --c 2.99792458e10 --v 2.99792458e8 --rho 1.2 --kappa 83333.33333333333 \
        --drift-tol 1 --quiet
  done
  for n in t6_n64 t6_n128; do
    d=$R/${F}_$n; printf "%-5s %-12s " $F $n
    python3 t6_marshak.py $d/tab/*.m1.00001.tab --hydro $d/tab/*.hydro_w.00001.tab $T6
  done
done
