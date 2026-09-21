#!/bin/bash -l
# Re-derives every number in the milestone-1c gate tables from the run directories.
cd "$(dirname "$0")/.."
A5="--c 2.99792458e10 --rho 1.2 --kappa 500 --v 3.0e7 --drift-tol 1 --quiet"
A1="--c 2.99792458e10 --rho 1.2 --kappa 100 --v 1.0e6 --drift-tol 1 --quiet"
R=runs_1c
echo "### G2/item1  T3 (static thick pulse), ap_form = unified vs alpha2"
for F in unified alpha2; do
  for n in tau1e1_dc tau1e1_plm tau1e3_dc tau1e3_plm tau1e6_dc tau1e6_plm \
           tau1e1_plm_n256; do
    d=$R/t3_${F}_${n}
    k=$(grep -oP 'kappa_s=\K[0-9.e+]+' $d/run.log|head -1)
    printf "%-8s %-18s " $F $n
    python3 t3_pulse.py $d/bin/*.bin --c 1 --rho 1 --kappa "$k" --subtract-min --quiet
  done
  printf "%-8s %-18s " $F nyquist
  python3 t3_pulse.py $R/t3_${F}_nyq_plm/bin/*.bin --c 1 --rho 1 --kappa 128000 \
      --nyquist --quiet
  printf "%-8s %-18s " $F tophat
  python3 t3_pulse.py $R/t3_${F}_tophat/bin/*.bin --c 1 --rho 1 --kappa 1280000 \
      --subtract-min --quiet
done
echo
echo "### G6  T3b redefined, case A: tau_cell 1 -> 1e3 (t = 1e5)"
for n in ap_unified ap_alpha2 none scaled; do
  printf "%-12s " $n
  python3 t3b_jump.py $R/t3bA_$n/tab/*.00002.tab --c 1 --kappa 64.0 --flux 1.0e-6 --quiet
done
echo "### G6  T3b redefined, case B: tau_cell 1e-3 -> 1 (t = 2000)"
for n in B_ap_unified B_ap_alpha2 B_none B_scaled; do
  printf "%-14s " $n
  python3 t3b_jump.py $R/t3b_$n/tab/*.00002.tab --c 1 --kappa 0.064 --flux 1.0e-3 --quiet
done
echo
echo "### G3  T4, dynamic case (kappa 500, v = 3e7, beta tau = 14)"
S5=$R/t4_k500_v0_512/tab/m1_advect_pulse.m1.00004.tab
S5k=$R/t4_k500_v0_1024/tab/m1_advect_pulse.m1.00004.tab
SA=$R/t4_a2_k500_v0_512/tab/m1_advect_pulse.m1.00004.tab
for n in d512_aphll_split d512_aphll_nosplit d512_scaled_split d512_none_split \
         d512_none_nosplit d512_aphll_split_ovc d512_aphll_split_nsub1; do
  d=$R/t4_$n; printf "%-24s " $n
  python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab $A5 --static $S5
done
for n in a2_d512_split a2_d512_nosplit; do
  d=$R/t4_$n; printf "%-24s " $n
  python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab $A5 --static $SA
done
for n in d1024_aphll_split d1024_aphll_nosplit; do
  d=$R/t4_$n; printf "%-24s " $n
  python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab $A5 --static $S5k
done
echo "### G3  T4, static case (kappa 100, v = 1e6, beta = 3.3e-5)"
S1=$R/t4_s512_aphll_split/tab/m1_advect_pulse.m1.00004.tab
S1k=$R/t4_k100_v0_1024b/tab/m1_advect_pulse.m1.00004.tab
d=$R/t4_k100_v1e6_512;  printf "%-24s " k100_v1e6_512
python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab $A1 --static $S1
d=$R/t4_k100_v1e6_1024; printf "%-24s " k100_v1e6_1024
python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab $A1 --static $S1k
echo
echo "### G4  T4b, uniform medium, beta tau_cell = 1e3"
for n in aphll_split_full aphll_split_ovc aphll_nosplit_full aphll_nosplit_ovc \
         none_split_full scaled_split_full; do
  d=$R/t4b_$n; printf "%-20s " $n
  python3 t4_advect.py $d/tab/*.m1.*.tab --hydro $d/tab/*.hydro_w.*.tab \
      --c 2.99792458e10 --v 2.99792458e8 --rho 1.2 --kappa 83333.33333333333 \
      --drift-tol 1 --quiet
done
echo
echo "### G5  T6, constant-c_v Marshak wave, reference = runs_1c/t6_ref_sn.txt"
T6="--table $R/t6_ref_sn.txt --marshak --a-rad 1e30 --cv 1.5 --rho 1 --centre 0.0 --quiet"
for N in "" n64_ n32_; do
  for n in ap_unified ap_alpha2 none scaled; do
    d=$R/t6${N:+${N%_}}_$n; [ -d "$d" ] || d=$R/t6_$n
    printf "%-10s %-12s " "${N:-n128_}" $n
    python3 t6_marshak.py $d/tab/*.m1.00001.tab --hydro $d/tab/*.hydro_w.00001.tab $T6
  done
done
