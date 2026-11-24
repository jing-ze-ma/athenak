#!/bin/bash
# orion NOTE 23zg: m35mms dt diagnostic, restart (read in place) from prod_ia11_1123/m35mms/rst/m35mms.00010.rst (cycle 21945, t ~1.728e5), IA11 binary, hpc-sm-01-17, 4 ranks, 100 cycles
J=/resnick/groups/carnegie_poc/jingze; D=$J/diag_m35dt_1123; BIN=$J/wt_orion1121/bin/athena_ia11_bd9553537; RST=$J/prod_ia11_1123/m35mms/rst/m35mms.00010.rst; MW=$J/orion_runs_1121/ctl/memwrap.sh
P=$J/prod_fresh_1121b; while IFS= read -r l; do export "$l" 2>/dev/null; done < $P/i1_m35mms_f1k.env; export VGD_XFE=0 ATHENA_DT_REPORT=10
cd $D; echo "$(date '+%F %T') start BIN md5=$(cut -c1-32 $BIN.md5) RST=$RST size=$(stat -c %s $RST)" >> $D/run.log
OMPI_MCA_io=romio321 srun --overlap --jobid=7333072 -w hpc-sm-01-17 -N 1 -n 4 --mpi=pmix -t 40 -J d_m35dt \
  env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 ATHENA_DT_REPORT=10 $MW $D $BIN -r $RST -d $D \
  time/nlim=22045 job/final_restart=false output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 < /dev/null > $D/log.0 2>&1
echo "$(date '+%F %T') rc=$? last: $(grep -E '^elapsed=.*cycle=' $D/log.0 | tail -1)" >> $D/run.log
