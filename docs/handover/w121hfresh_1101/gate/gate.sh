#!/bin/bash -l
# w121hfresh_1101 1x CPU gates (login node, nice, 12 ranks): parse + 2-cycle fresh start with ghosts.
# CPU binary athena_dhj_cpu_8bcd50c8 = fork/fin51-r-fix 182477553 + dhj-remap-mhd merge (as toi1518b gate).
module purge; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
export OMP_NUM_THREADS=1
D=/viper/ptmp2/jinma/w121hfresh_1101; G=$D/gate1x
B=/viper/ptmp2/jinma/builds/bin/athena_dhj_cpu_8bcd50c8
DECK=$D/deck/hfresh_1x.athinput
cd $G || exit 1
md5sum $B $DECK
$B -n -i $DECK > $G/parse.log 2>&1; echo "parse rc=$?"
sed 's/^variable *= hydro_w .*/&\nghost_zones = true/' $DECK > $G/deck_gz.athinput
K2="time/ndiag=1 output1/dt=1.0e-6 output2/dt=1.0e-6 output3/dt=1.0e9 output4/dt=1.0e9"
d=$G/run_fresh; mkdir -p $d; cd $d || exit 1
ln -sfn $D/pkg_1x/ck ck; ln -sfn $D/pkg_1x/ic ic
nice mpirun -np 12 --oversubscribe $B -i $G/deck_gz.athinput $K2 time/nlim=2 > run.log 2>&1; echo "fresh rc=$?"
echo GATE_DONE
