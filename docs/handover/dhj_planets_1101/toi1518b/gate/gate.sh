#!/bin/bash -l
# TOI-1518b CPU gates (login node, nice, 12 MPI ranks):
#   1. parse (-n) of the deck;
#   2. fresh start of the production deck (C32 x nx1 120, 24 MeshBlocks), ghost zones written,
#      NLIM cycles (default 2), counters every cycle, t = 0 bin -> gate_ic.py.
# CPU binary: athena_dhj_cpu_8bcd50c8 = fork/fin51-r-fix 182477553 + the dhj-remap-mhd merge
# (MHD remap pgen path, default off; rad_m1 halo defaults: no rad_m1 here).
module purge; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
export OMP_NUM_THREADS=1
P=/viper/ptmp2/jinma/dhj_planets_1101/toi1518b; G=$P/gate
B=/viper/ptmp2/jinma/builds/bin/athena_dhj_cpu_8bcd50c8
DECK=${DECK:-$P/deck/toi1518b_hfresh.athinput}
NLIM=${NLIM:-2}
cd $G || exit 1
$B -n -i $DECK > $G/${PARSE:-parse.log} 2>&1; echo "parse rc=$?"
sed 's/^variable *= hydro_w .*/&\nghost_zones = true/' $DECK > $G/${GZ:-deck_gz.athinput}
K2="time/ndiag=1 output1/dt=1.0e-6 output2/dt=1.0e-6 output3/dt=1.0e9 output4/dt=1.0e9"
d=$G/${RUN:-run_fresh}; mkdir -p $d; cd $d || exit 1
ln -sfn $P/ck ck; ln -sfn $P/ic ic
nice mpirun -np 12 --oversubscribe $B -i $G/${GZ:-deck_gz.athinput} $K2 time/nlim=$NLIM \
  > run.log 2>&1; echo "fresh rc=$?"
echo GATE_DONE
