#!/bin/bash -l
# steady states, be named vs hesdirk2 (the new default, key absent), run long enough
W=/viper/ptmp2/jinma/sph2_0924; I=$W/inp; X=$W/bin/athena_v2_none_cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
run() { local n=$1; shift; local d=$W/cpu/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np 1 $X -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) & }
T="rad_m1/implicit_tol=1.0e-10 rad_m1/implicit_lin_tol=1.0e-11"
for a in be h2; do s=$([ $a = be ] && echo _be || echo "")
  run s2fs_$a -i $I/sp_sph_fs$s.athinput $T time/nlim=200 output1/dcycle=100
  run s2atm_$a -i $I/sp_sph_atm$s.athinput mesh/nx1=32 meshblock/nx1=32 time/nlim=2000 output1/dcycle=500
  run s2vatm_$a -i $I/sp_sph_atm_vc$s.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=3000 output1/dcycle=500
  run s2vatmq_$a -i $I/sp_sph_atm_vc$s.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=3000 output1/dcycle=500 rad_m1/vet_col_surface_q=true
  run s2vstr_$a -i $I/sp_sph_atm_str_vc$s.athinput mesh/nx1=32 meshblock/nx1=32 time/nlim=3000 output1/dcycle=500
  run s2vpp_$a -i $I/milne_vc$s.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=3000 output1/dcycle=500
done
wait
for c in fs atm vatm vatmq vstr vpp; do
  fs=($(ls $W/cpu/s2${c}_h2/tab/*m1*.tab)); n=${#fs[@]}
  echo "SS2 $c h2 vs be: $(python3 $W/scripts/ana.py steady 'tab/*m1*.tab' $W/cpu/s2${c}_h2 $W/cpu/s2${c}_be)"
  echo "    be last two dumps: $(python3 $W/scripts/ana.py steadyf ${fs[$((n-2))]/_h2/_be} ${fs[$((n-1))]/_h2/_be})"
  echo "    $(grep -ho 'NON-CONVERGED=[^ ]*' $W/cpu/s2${c}_be/log.txt $W/cpu/s2${c}_h2/log.txt | tr '\n' ' ') $(grep -h 'time_scheme=hesdirk2:' $W/cpu/s2${c}_h2/log.txt)"
done
echo SS2 DONE
