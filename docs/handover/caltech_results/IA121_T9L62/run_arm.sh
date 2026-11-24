#!/bin/bash
# run_arm.sh S|L: orion NOTE 23zn t9l62 arms, athena_t9arm_932dbd534, from prod_ia11_1123/t9l62/rst_keep/t9l62.00012.rst (read only), N2 7333073 node B hpc-sm-02-17, 4 ranks (8 blocks, 2/GPU)
# keys = IA121_t9l62_restart.args exactly with ONE change (S: hydro/floor_pay_skip=false; L: hydro/dfloor_flux_limit=false) + nlim 6850, no rst/bin output
J=/resnick/groups/carnegie_poc/jingze; V=$1; D=$J/t9arms_1123/$V; BIN=$J/wt_orion1121/bin/athena_t9arm_932dbd534; RST=$J/prod_ia11_1123/t9l62/rst_keep/t9l62.00012.rst; MW=$J/orion_runs_1121/ctl/memwrap.sh
ARGS=$J/prod_ia121_1123/IA121_t9l62_restart.args; KEYS=$(grep -v '^#' $ARGS | tr '\n' ' ')
case $V in S) KEYS=${KEYS/hydro\/floor_pay_skip=true/hydro\/floor_pay_skip=false};; L) KEYS=${KEYS/hydro\/dfloor_flux_limit=true/hydro\/dfloor_flux_limit=false};; *) exit 1;; esac
KEYS="$KEYS time/nlim=6850 job/final_restart=false output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30"
P=$J/prod_fresh_1121b; while IFS= read -r l; do export "$l" 2>/dev/null; done < $P/i1_m35mms_f1k.env; export VGD_XFE=0
mkdir -p $D; cd $D; echo "$(date '+%F %T') start $V BIN md5=$(cut -c1-32 $BIN.md5) args md5=$(md5sum < $ARGS | cut -c1-32) RST=$RST KEYS: $KEYS" >> $D/run.log
OMPI_MCA_io=romio321 srun --overlap --jobid=7333073 -w hpc-sm-02-17 -N 1 -n 4 --mpi=pmix -t 180 -J t9arm_$V \
  env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 $MW $D $BIN -r $RST -d $D $KEYS < /dev/null > $D/log.0 2>&1
echo "$(date '+%F %T') rc=$? last: $(grep -E '^elapsed=.*cycle=' $D/log.0 | tail -1)" >> $D/run.log
