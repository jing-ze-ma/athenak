#!/bin/bash -l
# vet_col under hesdirk2: D^n (default) vs time2_vet_extrap = true, the TV2 relaxation
W=/viper/ptmp2/jinma/sph2_0924; I=$W/inp; X=$W/bin/athena_v2_none_cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
O="time/nlim=-1 time/tlim=6.4 rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=0.0 rad_m1/kappa_s=100.0 rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_lin_tol=1.0e-13"
run() { local n=$1; shift; local d=$W/cpu/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np 1 $X -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) & }
for c in 0.3 0.15 0.075 0.0375 0.009375; do
  run twx_$c -i $I/atmvc_h2x.athinput $O time/cfl_number=$c
  run twxq_$c -i $I/atmvc_h2x.athinput $O time/cfl_number=$c rad_m1/vet_col_surface_q=true
  run twq_$c -i $I/atmvc_h2.athinput $O time/cfl_number=$c rad_m1/vet_col_surface_q=true
done
wait
for a in twx twxq twq; do echo "TVX $a: $(python3 $W/scripts/ana.py order 'tab/sa.m1.*.tab' m1_e $W/cpu/${a}_0.009375 $W/cpu/${a}_{0.3,0.15,0.075,0.0375})"; done
echo "TVX twx vs tw_ref (extrap vs D^n at the finest dt): $(python3 $W/scripts/ana.py steady 'tab/sa.m1.*.tab' $W/cpu/twx_0.009375 $W/cpu/tw_ref)"
grep -h "NON-CONVERGED\|time_scheme=hesdirk2:" $W/cpu/twx_0.3/log.txt
echo TVX DONE
