#!/bin/bash
# usage: gate_sum.sh <ref> <arm> [<ref> <arm> ...]: rc, NON-CONVERGED, hst compare, bin files
W=/viper/ptmp2/jinma/m1int_0924/cpu
while [ $# -ge 2 ]; do
  a=$1; b=$2; shift 2
  nc=$(grep -o "NON-CONVERGED=[0-9.e+-]*" $W/$b/log.txt | tr "\n" " "); rc=$(tail -1 $W/$b/log.txt)
  nb=0; nd=0
  for f in $(cd $W/$a && find ./bin -name "*.bin" | sort | tail -4); do
    nb=$((nb+1)); python3 $W/../scripts/bincmp.py $W/$a/$f $W/$b/$f || nd=$((nd+1)); done
  echo "$b vs $a: $rc $nc lastbins_differ=$nd/$nb :: $(python3 $W/../scripts/cmp.py $W $a $b | cut -d' ' -f2-)"
done
