#!/bin/bash -l
# T-S4 (runs_5d/5e grey spherical atmosphere, vet_col + surface_q, be, steady):
# old = ref binary, new = <tag> binary with implicit_marshak_face = linear and
# vet_col_order2 = true.  usage: ts4.sh <tag> [ncore]
W=/viper/ptmp2/jinma/sporder2_0925; t=$1; nc=${2:-8}
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
for n in 32 64 128 256; do
  L=800; [ $n = 256 ] && L=1600
  for a in old new; do
    b=$W/bin/athena_${t}_none_cpu; [ $a = old ] && b=$W/bin/athena_ref_none_cpu
    d=$W/cpu/ts4/${a}${n}_nc$nc; rm -rf $d; mkdir -p $d
    (cd $d && OMP_NUM_THREADS=1 nice mpirun -np 1 $b -i $W/inp/ts4_$a.athinput -d $d \
      mesh/nx1=$n meshblock/nx1=$n time/nlim=$L rad_m1/vet_col_ncore=$nc > log.txt 2>&1) &
  done
done
wait
echo TS4 DONE
