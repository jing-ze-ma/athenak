#!/bin/bash -l
# lin_tol accuracy on the CPU: He slab (vet_sc, 100 cycles, 1 rank) and He box
# (84x32x32 vet_sc, 30 cycles, 4 ranks) at implicit_lin_tol 1e-10 / 1e-9 / 1e-8 against
# a TIGHT reference (implicit_tol 1e-12, lin_tol 1e-11).  usage: ltol_cpu.sh <exe> <dir>
E=$1; R=$2; G=/viper/ptmp2/jinma/wt_fast3/tests_m1/gates
source /etc/profile.d/modules.sh 2>/dev/null; module purge >/dev/null 2>&1
module load gcc/14 openmpi/5.0 >/dev/null 2>&1
BOX="box3d_be_x mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=30 rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
SLAB="slab2d_plm_vimp_be_x time/nlim=100 rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
TIGHT="rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_lin_tol=1.0e-11"
go() {  # go <name> <np> <geo...>
  local n=$1 np=$2; shift 2; local inp=$1; shift
  local d=$R/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none $E \
     -i $G/$inp.athinput -d $d problem/vpert=1.0e-2 "$@" > log.txt 2>&1) &
}
for g in box slab; do
  np=4; geo=$BOX; [ $g = slab ] && { np=1; geo=$SLAB; }
  go ${g}_ref $np $geo $TIGHT
  go ${g}_l10 $np $geo rad_m1/implicit_lin_tol=1.0e-10
  go ${g}_l9 $np $geo rad_m1/implicit_lin_tol=1.0e-9
  go ${g}_l8 $np $geo rad_m1/implicit_lin_tol=1.0e-8
done
wait
cd $G
for g in box slab; do for a in l10 l9 l8; do
  python3 -c "
import cmp,sys
d=cmp.diff('$R/${g}_$a','$R/${g}_ref')
print('$g $a', ' '.join('%s=%s'%(k,v) for k,v in d.items()))"
  grep -o "NON-CONVERGED=[^ ]*" $R/${g}_$a/log.txt | tail -1
  grep -o "inner iterations mean=[^ ]*" $R/${g}_$a/log.txt | tail -1
done; done
