#!/bin/bash -l
# held-gas VET steady states: be vs hesdirk2 (base = time2_vet_gas old, new = start)
W=/viper/ptmp2/jinma/h2acc_0924; Q=/viper/ptmp2/jinma/sph2_0924/inp; S=$W/scripts
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
run() { local n=$1 x=$2; shift 2; local d=$W/cpu/ss/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np 1 $W/bin/athena_${x}_none_cpu -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) & }
for cf in 0.3 0.075; do
  f=$(python3 -c "print(int(round(0.3/$cf)))")
  for a in be h2b h2n; do
    case $a in be) s=_be; x=base;; h2b) s=""; x=base;; h2n) s=""; x=new;; esac
    C="time/cfl_number=$cf"
    run vatm_${cf}_$a $x -i $Q/sp_sph_atm_vc$s.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=$((3000*f)) output1/dcycle=$((500*f)) $C
    run vstr_${cf}_$a $x -i $Q/sp_sph_atm_str_vc$s.athinput mesh/nx1=32 meshblock/nx1=32 time/nlim=$((3000*f)) output1/dcycle=$((500*f)) $C
    run vpp_${cf}_$a $x -i $Q/milne_vc$s.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=$((3000*f)) output1/dcycle=$((500*f)) $C
    run vsc_${cf}_$a $x -i $Q/milne_vc$s.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=$((3000*f)) output1/dcycle=$((500*f)) $C rad_m1/closure=vet_sc
  done
done
wait
for cf in 0.3 0.075; do for c in vatm vstr vpp vsc; do
  echo "SS $c cfl=$cf base-h2 vs be: $(python3 $S/ana.py steady 'tab/*m1*.tab' $W/cpu/ss/${c}_${cf}_h2b $W/cpu/ss/${c}_${cf}_be 2>/dev/null | tail -n1)"
  echo "SS $c cfl=$cf new-h2  vs be: $(python3 $S/ana.py steady 'tab/*m1*.tab' $W/cpu/ss/${c}_${cf}_h2n $W/cpu/ss/${c}_${cf}_be 2>/dev/null | tail -n1)"
  fs=($(ls $W/cpu/ss/${c}_${cf}_be/tab/*m1*.tab)); n=${#fs[@]}
  echo "    be steadiness: $(python3 $S/ana.py steadyf ${fs[$((n-2))]} ${fs[$((n-1))]} 2>/dev/null)"
  echo "    $(grep -ho 'NON-CONVERGED=[^ ]*' $W/cpu/ss/${c}_${cf}_{be,h2b,h2n}/log.txt | tr '\n' ' ') $(grep -ho 'stage fallbacks=[^ ]*' $W/cpu/ss/${c}_${cf}_h2n/log.txt) $(tail -qn1 $W/cpu/ss/${c}_${cf}_{be,h2b,h2n}/log.txt | tr '\n' ' ')"
done; done
echo SS DONE
