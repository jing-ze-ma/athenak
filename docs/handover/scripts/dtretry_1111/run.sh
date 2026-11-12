#!/bin/bash
# dtretry-1111 CPU gate (Orion): HG small cs deck (6 blocks 704x16x16) + HG production keys + seed, 4 cycles,
# restart dump at cycle 4; arms run concurrently, 3 MPI ranks each (21 cores).
# usage: BASE=<athena of 84cd20888> FIX=<athena of 3ec401877> TB=<dir with the 3 tables> ./run.sh
S=$(cd $(dirname $0) && pwd); R=${RUNS:-$PWD/runs}; mkdir -p $R
: ${BASE:?} ${FIX:?} ${TB:?}
export OMP_NUM_THREADS=1
PK=$(cat $S/prodkeys.txt)
SEED="problem/dbg_vac_seed=1e-9 problem/dbg_vac_rlo=6.52e12 problem/dbg_vac_rhi=7.15e12 problem/dbg_vac_frac=0.05 problem/dbg_vac_vr=1.5e8 problem/dbg_vac_omega=2e-7"
B="problem/he_ic_file=$TB/ic_hegiant_14p5M_mk.txt problem/he_opac_table=$TB/rosseland_mesakap_afe2dd8da2bd.txt problem/he_planck_table=$TB/planck_ext2_x0.0009496_z0.0215304_416a61196c64.txt output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 mesh/nx2=16 mesh/nx3=16"
T="rad_m1/implicit_hr_cmem=true rad_m1/implicit_hr_cmem_cap=1.0 rad_m1/time2_vet_x2_fields=15 rad_m1/implicit_physflux=true rad_m1/implicit_physflux_mom=fa rad_m1/implicit_physflux_zold=fc rad_m1/implicit_hr_cmem_rst=true rad_m1/implicit_hr_cmem_fainit=true"
F="time/dt_retry_force_fac=1.0 time/dt_retry_force_cycle=3"
P="time/dt_retry_views=all time/dt_retry_audit=true time/dt_retry_poison=1"
arm() { local d=$R/$1; mkdir -p $d; cd $d || return
  python3 $S/mkdeck.py $S/hg_small_deck.athinput $d/deck.athinput $PK $SEED output5/dcycle=4 $3
  timeout 7200 mpirun -n 3 --bind-to none $2 -i $d/deck.athinput $B time/nlim=4 -d $d < /dev/null > run.log 2>&1
  echo "$1 rc=$? fatal=$(grep -c FATAL run.log) postrig=$(grep -c POS_TRIGGER run.log) $(grep 'dt_retry poison:' run.log | cut -c1-60) $(grep -E '^cycle=' run.log | tail -1 | cut -c1-50)"; }
arm BN  $BASE ""      & arm CN  $FIX "" & arm CF  $FIX "$F" & arm CP  $FIX "$F $P" &
arm CTN $FIX "$T"     & arm CTF $FIX "$T $F" & arm CTP $FIX "$T $F $P" & wait
cd $R
for p in "BN CN" "CN CF" "CN CP" "CTN CTF" "CTN CTP"; do set -- $p; echo "== $1 vs $2"; python3 $S/cmp.py $1 $2; done
