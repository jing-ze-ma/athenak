#!/bin/bash -l
# vet_col steady state be vs hesdirk2: (a) pure scattering (the source is E, no gas T);
# (b) thermal at cfl 0.075 (does the 2e-4 of cfl 0.3 scale with dt?); (c) T-S2 at tol 1e-10
W=/viper/ptmp2/jinma/sph2_0924; I=$W/inp; X=$W/bin/athena_new_none_cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
run() { local n=$1; shift; local d=$W/cpu/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np 1 $X -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) & }
S="rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=0.0 rad_m1/kappa_s=100.0"
G="mesh/nx1=64 meshblock/nx1=64 output1/dcycle=500"
for a in be h2; do s=$([ $a = be ] && echo _be || echo "")
  run s3sc_$a -i $I/sp_sph_atm_vc$s.athinput $G time/nlim=3000 $S
  run s3c075_$a -i $I/sp_sph_atm_vc$s.athinput $G time/nlim=12000 output1/dcycle=2000 time/cfl_number=0.075
  run s3fs_$a -i $I/sp_sph_fs$s.athinput rad_m1/implicit_tol=1.0e-10
done
wait
for c in sc c075 fs; do
  echo "SS3 $c h2 vs be: $(python3 $W/scripts/ana.py steady 'tab/*.tab' $W/cpu/s3${c}_h2 $W/cpu/s3${c}_be) $(grep -ho 'NON-CONVERGED=[^ ]*' $W/cpu/s3${c}_be/log.txt $W/cpu/s3${c}_h2/log.txt | tr '\n' ' ') $(grep -h 'time_scheme=hesdirk2:' $W/cpu/s3${c}_h2/log.txt)"
done
echo SS3 DONE
