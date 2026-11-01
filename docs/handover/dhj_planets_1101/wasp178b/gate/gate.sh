#!/bin/bash
# CPU gates for one deck of wasp178b (login node, nice, 8 MPI ranks): -n parse; fresh start with ghost zones
# on a reduced angular mesh (nx2 = nx3 = 16, MeshBlocks 8x8: same radial grid, same IC fill), 2 cycles.
# usage: gate.sh DECK NAME [extra keys]   (run dirs gate/run_NAME, parse log gate/parse_NAME.log)
# binary: fin51-r-fix 182477553 CPU build (= the source of the hfresh GPU binary), built by the
#   kelt20b worker, /viper/ptmp2/jinma/dhj_planets_1101/kelt20b_work/athena_dhj_cpu_18247755 (used read-only)
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false
P=/viper/ptmp2/jinma/dhj_planets_1101/wasp178b; G=$P/gate; N=/viper/ptmp2/jinma/dhj_planets_1101/kelt20b_work/athena_dhj_cpu_18247755
DECK=$1; NAME=$2; shift 2
md5sum $N $DECK
d=$G/run_$NAME; mkdir -p $d; cd $d || exit 1
ln -sfn $P/ck ck; ln -sfn $P/ic ic
nice $N -n -i $DECK "$@" > $G/parse_$NAME.log 2>&1; echo "parse $NAME rc=$?"
K2="mesh/nx2=16 mesh/nx3=16 meshblock/nx2=8 meshblock/nx3=8 time/ndiag=1 output1/dt=1.0e-6 output2/dt=1.0e-6 output3/dt=1.0e9 output4/dt=1.0e9"
sed 's/^variable *= hydro_w .*/&\nghost_zones = true/' $DECK > $G/deck_gz_$NAME.athinput
nice mpirun -np 8 --bind-to none $N -i $G/deck_gz_$NAME.athinput $K2 time/nlim=2 "$@" > run.log 2>&1; echo "fresh $NAME rc=$?"
grep -c FATAL run.log $G/parse_$NAME.log
echo GATE_DONE
