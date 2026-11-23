#!/bin/bash
# run_diag.sh HOLDER NODE A|B: orion NOTE 23q/23s IA9 t9l62 diagnostic restart (IA9 binary, kept rst t9l62.00001 read in place, 300 cycles 6450 -> 6750), 1 node 4 ranks
# A = 23q keys; B = A + rad_m1/time2_k1_xpat=pos (23s). Run dir $D/<V>.
J=/resnick/groups/carnegie_poc/jingze; D0=$J/ia9_t9diag_1123; HJ=$1; ND=$2; V=$3; : "${HJ:?}" "${ND:?}" "${V:?}"; D=$D0/$V; mkdir -p $D
case $V in A) EX="";; B) EX="rad_m1/time2_k1_xpat=pos";; *) echo "bad variant"; exit 1;; esac
BIN=$J/wt_orion1121/bin/athena_ia9_77551e823; RST=$J/prod_ia9_1123/t9l62/rst/t9l62.00001.rst; MW=$J/orion_runs_1121/ctl/memwrap.sh
P=$J/prod_fresh_1121b; while IFS= read -r l; do export "$l" 2>/dev/null; done < $P/i1_m35mms_f1k.env; export VGD_XFE=0
cd $D; echo "$(date '+%F %T') start $V HOLDER=$HJ NODE=$ND BIN md5=$(cut -c1-32 $BIN.md5) RST=$RST EX=hydro/dt_report=1 hydro/floor_audit=1 time/nlim=6750 $EX" >> $D/run.log
OMPI_MCA_io=romio321 srun --overlap --jobid=$HJ -w $ND -N 1 -n 4 --mpi=pmix -t 60 -J d_t9$V \
  env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 $MW $D $BIN -r $RST -i $D0/deck_diag.athinput -d $D \
  hydro/dt_report=1 hydro/floor_audit=1 time/nlim=6750 $EX < /dev/null > $D/log.0 2>&1
echo "$(date '+%F %T') rc=$? last: $(grep -E '^elapsed=.*cycle=' $D/log.0 | tail -1)" >> $D/run.log
