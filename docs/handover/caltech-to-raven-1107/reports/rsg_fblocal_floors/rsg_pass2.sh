#!/bin/bash -l
# rsg2.sh: RSG real state, same rst copy as rsg.sh, on vs off, cs_cross_mono + guard_mode 7 + implicit_fb_local + implicit_fb_diag, until END (default 04:45 PST), arms swapped relative to pass 1 so node bias shows
J=/resnick/groups/carnegie_poc/jingze; E=$J/cstan_1106/real; P=$J/rsg_E_1020; JOB=5691175; NB=${1:-5b03746b2}; ABIN=$J/builds_rstfixx/bin/athena_gpu_he_star_m1_cts_${NB}_fmadoff
XK=${XK:-"rad_m1/cs_cross_mono=true rad_m1/implicit_opac_newton_guard_mode=7 rad_m1/implicit_fb_local=true rad_m1/implicit_fb_diag=200"}
source $P/runH.cfg; BIN=$ABIN; RS=$E/rst/rsg.rst; INFO=$J/bsg_1009/files/rst_info.py
END=${END:-"2026-11-07 04:45:00"}; S=$(( $(date -d "$END" +%s) - $(date +%s) - 300 )); [ $S -lt 600 ] && { echo "too late"; exit 1; }; WL=$(printf '%02d:%02d:%02d' $((S/3600)) $(((S%3600)/60)) $((S%60)))
INF=$(srun --overlap --jobid=$JOB -w hpc-sm-02-09 -N1 -n1 -c1 --mem=8G python3 $INFO $RS); OK=$(echo "$INF" | tail -1); CY=$(srun --overlap --jobid=$JOB -w hpc-sm-02-09 -N1 -n1 -c1 --mem=8G python3 $J/layinv_1105/gpu/rstcyc.py $RS | cut -d' ' -f1)
OUTK="output2/bin_double=true output2/dcycle=500 output3/bin_double=true output3/dcycle=500 output5/dcycle=1000000"
arm(){ local D=$1 NODE=$2 V=$3; rm -rf $D; mkdir -p $D; cd $D; echo "start $(date -Iseconds) $NODE tan=$V cycle0=$CY md5=$(md5sum < $BIN | cut -c1-32) wall=$WL keys=$XK" > arm.info
  timeout $((S+600)) srun --overlap --jobid=$JOB -w $NODE -N1 -n 4 -c 8 --cpu-bind=none --mpi=pmix --kill-on-bad-exit=1 env VGD_TANQ_DIAG=20 $E/wrap.sh $BIN -r $RS $OK $XKEYS $OUTK $XK rad_m1/vet_gd_tangential=$V time/nlim=$((CY+100000)) -d $D -t $WL < /dev/null > run.log 2>&1; echo "end $(date -Iseconds) rc=$?" >> arm.info; }
arm $E/rsg_on3 hpc-sm-02-09 true & arm $E/rsg_off3 hpc-sm-02-16 false & wait
