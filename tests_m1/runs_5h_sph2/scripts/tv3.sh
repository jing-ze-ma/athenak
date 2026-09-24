#!/bin/bash -l
# TV2 setup (pure scattering relaxation of the T-S4 atmosphere, t = 6.4) with eddington
# and m1 under hesdirk2: does the fixed-D^n vet_col tensor limit the order?
W=/viper/ptmp2/jinma/sph2_0924; I=$W/inp; X=$W/bin/athena_new_none_cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
O="time/nlim=-1 time/tlim=6.4 rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=0.0 rad_m1/kappa_s=100.0 rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_lin_tol=1.0e-13"
run() { local n=$1; shift; local d=$W/cpu/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np 1 $X -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) & }
for cl in eddington m1; do
  for c in 0.3 0.15 0.075 0.0375 0.009375; do run tw${cl}_h2_$c -i $I/atmvc_h2.athinput $O time/cfl_number=$c rad_m1/closure=$cl; done
done
for c in 0.3 0.15 0.075 0.0375; do run twtl_h2_$c -i $I/atmvc_h2.athinput $O time/cfl_number=$c rad_m1/implicit_tol=1.0e-14 rad_m1/implicit_lin_tol=1.0e-15; done
run twtl_ref -i $I/atmvc_h2.athinput $O time/cfl_number=0.009375 rad_m1/implicit_tol=1.0e-14 rad_m1/implicit_lin_tol=1.0e-15
wait
for cl in eddington m1; do echo "TV3 $cl h2: $(python3 $W/scripts/ana.py order 'tab/sa.m1.*.tab' m1_e $W/cpu/tw${cl}_h2_0.009375 $W/cpu/tw${cl}_h2_{0.3,0.15,0.075,0.0375})"; done
echo "TV3 vet_col h2 tol 1e-14: $(python3 $W/scripts/ana.py order 'tab/sa.m1.*.tab' m1_e $W/cpu/twtl_ref $W/cpu/twtl_h2_{0.3,0.15,0.075,0.0375})"
grep -h "NON-CONVERGED" $W/cpu/twtl_h2_0.3/log.txt
echo TV3 DONE
