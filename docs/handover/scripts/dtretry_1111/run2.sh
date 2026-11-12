#!/bin/bash
# dtretry-1111 CPU gate (Orion): HG small cs deck (6 blocks 704x16x16) + HG production keys + seed, 4 cycles,
# restart dump at cycle 4; arms run concurrently, 3 MPI ranks each (21 cores).
# usage: FIX=<athena of 34d198efa> TB=<dir with the 3 tables> ./run2.sh   (10 arms x 3 ranks = 30 cores)
S=$(cd $(dirname $0) && pwd); R=${RUNS:-$PWD/runs}; mkdir -p $R
: ${FIX:?} ${TB:?}
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
  echo "$1 rc=$? fatal=$(grep -c FATAL run.log) $(grep -E '^cycle=' run.log | tail -1 | cut -c1-50)"; }
# addendum (11-12): FIX = 34d198efa (m1.vgd_hr in the minimal list); V* = CF + one candidate group restored
A="time/dt_retry_views_add="
G1=hyd.u1,hyd.uflx.x1f,hyd.uflx.x2f,hyd.uflx.x3f,hyd.fofc_cnt,hyd.utest,hyd.wbq0,hyd.sbf_s,hyd.wder,hyd.wtemp
G2=m1.pos_cnt_d,m1.cf_r,m1.thin_frz,m1.ktn,m1.ost,m1.iw,m1.ifw,m1.ifw2,m1.ifw3,m1.f0x1n,m1.f0x2n,m1.f0x3n
G3=m1.gcz1,m1.gcz2,m1.gcz3,m1.vet_now,m1.vet_opac,m1.thx2,m1.thx3,m1.u1,m1.vcol_q,m1.thw,m1.thq,m1.krw
G4=m1.t2gsv,m1.t2gsv_h,m1.t2k2,m1.t2inc,m1.t2f1,m1.t2f2,m1.t2f3,m1.dbrow,m1.vmw,m1.d3w,m1.amw
G5=m1.vlat_lx,m1.vra_a0,m1.vra_a1,m1.vra_c,m1.vgd_rg_store,m1.vgd_cs,m1.vgd_cssb,m1.vgd_csrb,m1.vgd_csb,m1.vgd_crb
G6=m1.hm_sbuf,m1.hm_rbuf,m1.hm_shost,m1.hm_rhost,m1.lay_part,m1.lay_part_h,m1.vgd_sbuf,m1.vgd_rbuf,m1.vgd_bbuf0,am.amp_w,am.amp_u,am.amq_d
arm CN $FIX "" & arm CF $FIX "$F" & arm CTN $FIX "$T" & arm CTF $FIX "$T $F" &
for g in 1 2 3 4 5 6; do eval "arm V$g \$FIX \"\$F $A\$G$g\"" & done; wait
cd $R
for p in "CN CF" "CTN CTF" "CN V1" "CN V2" "CN V3" "CN V4" "CN V5" "CN V6"; do set -- $p; echo "== $1 vs $2"; python3 $S/cmp.py $1 $2; done
