#!/bin/bash
# G0 for the pgen switch: rad_m1_beam.athinput (beam_packet not named) with the d0c59f7c
# binary (build_mpi) and the m1-validate binary (build_mpi2); all dumps must cmp equal.
V=/viper/ptmp2/jinma/validate_0923; B=$V/bitwise
until [ -f $V/wt/build_mpi2.done ]; do sleep 20; done
cat $V/wt/build_mpi2.done
for b in build_mpi build_mpi2; do d=$B/$b; rm -rf $d; mkdir -p $d; cd $d
  nice -n 10 $V/wt/$b/src/athena -i $V/wt/inputs/tests/rad_m1_beam.athinput -d $d time/tlim=0.5 > log 2>&1; done
n=0; bad=0; for f in $B/build_mpi/bin/*.bin; do n=$((n+1)); cmp -s $f $B/build_mpi2/bin/$(basename $f) || bad=$((bad+1)); done
echo "bitwise: $n dumps compared, $bad differ"
