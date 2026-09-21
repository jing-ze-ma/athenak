#!/bin/bash -l
# Re-derives every number in the milestone-1c-B gate tables from the run directories.
# Run from anywhere; the dumps must still be present (they are deleted after
# RESULTS.txt and tests_m1/plots/*.json have been written).
cd "$(dirname "$0")/.."
R=runs_1cB

echo "### Step 0  the 1c-A defaults on the post-merge binary"
python3 t3_pulse.py $R/step0/t3_tau1e3_plm/bin/*.bin --c 1 --rho 1 --kappa 128000 \
    --subtract-min --quiet
A5="--c 2.99792458e10 --rho 1.2 --kappa 500 --v 3.0e7 --drift-tol 1 --quiet"
S=$R/step0/t4_s512/tab/m1_advect_pulse.m1.00004.tab
for n in split nosplit; do
  d=$R/step0/t4_d512_$n; printf "%-10s " $n
  python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab $A5 --static $S
done
for n in full ovc; do
  d=$R/step0/t4b_$n; printf "%-6s " $n
  python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab \
      --c 2.99792458e10 --v 2.99792458e8 --rho 1.2 --kappa 83333.33333333333 \
      --drift-tol 1 --quiet
done

echo
echo "### G-recon  pulse1d convergence (N = 32..512)"
for rec in plm ppm4 ppmx wenoz; do
  A=""
  for N in 32 64 128 256 512; do
    d=$R/p1d_${rec}_$N/tab
    A="$A --run $N $d/m1_pulse1d.m1.00000.tab $d/m1_pulse1d.m1.00001.tab"
  done
  printf "%-6s " $rec
  python3 t1_pulse.py $A --c 1.0 --quiet
done
echo "### G-recon  beam"
for rec in plm ppm4 ppmx wenoz; do
  printf "%-6s " $rec
  python3 t1_beam.py $R/beam_$rec/bin/m1_beam.m1.00006.bin --c 1 --quiet
done
echo "### G-recon  T3 rates"
for rec in plm ppm4 ppmx wenoz; do
  for t in tau1e1:1280 tau1e3:128000 tau1e6:128000000; do
    n=${t%%:*}; k=${t##*:}; printf "%-6s %-8s " $rec $n
    python3 t3_pulse.py $R/t3_${rec}_${n}/bin/*.bin --c 1 --rho 1 --kappa $k \
        --subtract-min --quiet
  done
  printf "%-6s %-8s " $rec nyquist
  python3 t3_pulse.py $R/t3_${rec}_nyq/bin/*.bin --c 1 --rho 1 --kappa 128000 \
      --nyquist --quiet
done
echo "### G-recon  T3 rate error vs N at tau_cell = 10"
for rec in plm ppm4 ppmx wenoz; do
  for N in 64 128 256 512; do
    K=$(python3 -c "print(1280.0*$N/128)")
    printf "%-6s N=%-4s " $rec $N
    python3 t3_pulse.py $R/t3n_${rec}_$N/bin/*.bin --c 1 --rho 1 --kappa $K \
        --subtract-min --quiet
  done
done

echo
echo "### G-shear  T4c, split_vel = recon vs cell, against the 2048-cell reference"
for v in recon cell; do
  printf "%-6s " $v
  python3 t4c_shear.py \
      --ref $R/shear_ref2048/tab/m1_advect_shear.m1.00001.tab \
      --run 64  $R/shear_${v}_64/tab/m1_advect_shear.m1.00001.tab \
      --run 128 $R/shear_${v}_128/tab/m1_advect_shear.m1.00001.tab \
      --run 256 $R/shear_${v}_256/tab/m1_advect_shear.m1.00001.tab \
      --run 512 $R/shear_${v}_512/tab/m1_advect_shear.m1.00001.tab --quiet --max-l1 1
done

echo
echo "### G-split-vel  T4 dynamic, split_vel = recon vs cell (the deciding gate)"
for N in 512 1024; do
  if [ $N = 512 ]; then SUF=""; ST=$R/step0/t4_s512; else SUF=k; ST=$R/recheck/t4k_static; fi
  for v in recon cell; do
    d=$R/recheck/t4${SUF}_$v; printf "%-5s %-6s " $N $v
    python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab $A5 \
        --static $ST/tab/m1_advect_pulse.m1.00004.tab
  done
done

echo
echo "### G-cons  T8c, 200 cycles of T4 dynamic"
d=$R/t8c
python3 t8_conserve.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab \
    --c 2.99792458e10 --tol 1e-12 --quiet
