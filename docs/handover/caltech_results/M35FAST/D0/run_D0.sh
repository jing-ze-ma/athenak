#!/bin/bash
# run_D0.sh [nomesh]: orion NOTE 23zm m35mms arm D0 (IA12.1 fa59bffd6) from prod_ia11_1123/m35mms/rst_keep/m35mms.00010.rst (read only), N2 7333073 node A hpc-sm-02-11, 4 ranks, 2000 cycles
J=/resnick/groups/carnegie_poc/jingze; D=$J/m35fast_1123/D0; BIN=$J/wt_orion1121/bin/athena_ia121_fa59bffd6; RST=$J/prod_ia11_1123/m35mms/rst_keep/m35mms.00010.rst; MW=$J/orion_runs_1121/ctl/memwrap.sh
KEYS=$(cat $J/m35fast_1123/D0_args.txt); [ "$1" = nomesh ] && KEYS=$(echo "$KEYS" | tr ' ' '\n' | grep -v '^mesh/' | tr '\n' ' ')
P=$J/prod_fresh_1121b; while IFS= read -r l; do export "$l" 2>/dev/null; done < $P/i1_m35mms_f1k.env; export VGD_XFE=0 ATHENA_DT_REPORT=25
cd $D; echo "$(date '+%F %T') start BIN md5=$(cut -c1-32 $BIN.md5) args md5=$(md5sum < $J/m35fast_1123/D0_args.txt | cut -c1-32) mode=${1:-full} RST=$RST" >> $D/run.log
OMPI_MCA_io=romio321 srun --overlap --jobid=7333073 -w hpc-sm-02-11 -N 1 -n 4 --mpi=pmix -J D0_m35 \
  env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 ATHENA_DT_REPORT=25 $MW $D $BIN -r $RST -d $D -t 05:00:00 $KEYS < /dev/null > $D/log.0 2>&1
echo "$(date '+%F %T') rc=$? last: $(grep -E '^elapsed=.*cycle=' $D/log.0 | tail -1)" >> $D/run.log
