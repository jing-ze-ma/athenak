#!/bin/bash -l
# rsg.sh [BIN tag]: RSG real state, copy of the newest rst, ON and OFF side by side (4 GPUs each, ~114 GB/GPU), arm A hpc-sm-02-09, arm B hpc-sm-02-16; then swapped. ONLY after the coordinator stopped RSG production.
J=/resnick/groups/carnegie_poc/jingze; E=$J/cstan_1106/real; P=$J/rsg_E_1020; JOB=5691175; NB=${1:-5b03746b2}; BIN=$J/builds_rstfixx/bin/athena_gpu_he_star_m1_cts_${NB}_fmadoff; NCYC=${NCYC:-3000}; WALL=${WALL:-02:40:00}; TMO=$((3600*3)); CCYC=${CCYC:-100}
XK=${XK:-"rad_m1/cs_cross_mono=true rad_m1/implicit_opac_newton_guard_mode=7"}
ABIN=$BIN; source $P/runH.cfg; BIN=$ABIN; mkdir -p $E/rst; RS=$E/rst/rsg.rst
if [ ! -f $RS ]; then N=$(ls $P/prodH/rst/*.rst | sort | tail -1); srun --overlap --jobid=$JOB -w hpc-sm-02-09 -N1 -n1 -c2 --mem=40G nice -n 19 cp $N $RS.tmp && mv $RS.tmp $RS; echo "$N" > $E/rst/rsg.src; fi
INF=$(srun --overlap --jobid=$JOB -w hpc-sm-02-09 -N1 -n1 -c1 --mem=8G python3 $INFO $RS); OK=$(echo "$INF" | tail -1); CY=$(srun --overlap --jobid=$JOB -w hpc-sm-02-09 -N1 -n1 -c1 --mem=8G python3 $J/layinv_1105/gpu/rstcyc.py $RS | cut -d' ' -f1)
OUTK="output2/bin_double=true output2/dcycle=500 output3/bin_double=true output3/dcycle=500 output5/dcycle=1000000"
arm(){ local NC=$NC WL=$WL D=$1 NODE=$2 V=$3; rm -rf $D; mkdir -p $D; cd $D; echo "start $(date -Iseconds) $NODE tan=$V cycle0=$CY md5=$(md5sum < $BIN | cut -c1-32) rst=$(cat $E/rst/rsg.src)" > arm.info
  timeout $TMO srun --overlap --jobid=$JOB -w $NODE -N1 -n 4 -c 8 --cpu-bind=none --mpi=pmix --kill-on-bad-exit=1 env VGD_TANQ_DIAG=20 $E/wrap.sh $BIN -r $RS $OK $XKEYS $OUTK $XK rad_m1/vet_gd_tangential=$V time/nlim=$((CY+NC)) -d $D -t $WL < /dev/null > run.log 2>&1; echo "end $(date -Iseconds) rc=$?" >> arm.info; }
NC=$NCYC; WL=$WALL; arm $E/rsg_off1 hpc-sm-02-09 false & arm $E/rsg_on1 hpc-sm-02-16 true & wait
NC=$CCYC; WL=00:20:00; TMO=1500; arm $E/rsg_on2 hpc-sm-02-09 true & arm $E/rsg_off2 hpc-sm-02-16 false & wait
