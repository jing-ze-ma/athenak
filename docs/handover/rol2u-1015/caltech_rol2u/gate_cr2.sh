#!/bin/bash -l
# gate_cr2.sh (CR2): FIN4-rst readability + BASE=PROD bitwise + ROL2U smoke/restart/forced-retry gate, from the SEED copy.
#   sbatch -J cr2g_$CASE -A carnegie_poc -p gpu --qos=debug -t 00:30:00 --nodes=1 --ntasks=$NGPU -c 8 \
#          --gres=gpu:nvidia_h200:$NGPU --mem=240G --exclude=hpc-sm-01-09 --export=ALL,CASE=hg|e1,BQ_MD5=<md5>,CR2_X=<bundle> gate_cr2.sh
# Arms (all -r SEED, nlim = seed cycle + NG; hst every cycle; production output cadence otherwise; final outputs written at the end):
#   P  production binary (d505d21f / 5b7e1bc5), production keys            -> the bitwise reference (cheap; production may be idle)
#   B  1bd0c2ce, production keys + time/dt_retry=true  (= BASE)            GATE-BASE: P = B BITWISE (rst data, final bins, hst rows)
#   C  1bd0c2ce + dt_retry + ROL2U (= the ROL2U chain's settings)
#   L1 = C to seed + NG/2;  L2 = -r L1 -> seed + NG                         GATE-RESTART: C = L2 BITWISE
#   F  = C + time/dt_retry_force_cycle=seed+5 time/dt_retry_force_fac=1.0   GATE-RETRY:   C = F  BITWISE (exactly 1 DT_RETRY line)
#   all: rc 0, FATAL 0, NaN 0; <par_end> offset of the C / L2 rst < 45056 (ParameterInput restart header limit) -> VERIFIED
set -u
source ${CR2_X:?export CR2_X=<this bundle dir>}/cr2_common.sh
chkmd5 $BQ $BQ_MD5; chkmd5 $BP $BP_MD5
S=$D/gate/j$SLURM_JOB_ID; mkdir -p $S; RST=$(seedrst); [ -f "$RST" ] || { echo "no seed rst: run seed_cr2.sh $CASE"; exit 1; }
read -r _ T0 NC0 _ <<< "$(seedinfo)"; K="$(seedkeys) $PKEYS output1/dt=1e-6"; N=$((NC0+NG)); NH=$((NC0+NG/2)); FC=$((NC0+5))
echo "CASE $CASE seed $RST t $T0 cycle $NC0 -> $N; node $(hostname); $(md5sum $BQ | cut -c1-8) $(md5sum $BP | cut -c1-8)"
run() { local n=$1 B=$2; shift 2; local R=$S/$n; mkdir -p $R; echo "$(md5sum $B | cut -d' ' -f1) $B $*" > $R/cmd.txt
  nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits -lms 500 > $R/mem.txt 2>/dev/null & local mp=$!
  local t0=$SECONDS; ( cd $R; $LAUNCH $B "$@" -d $R < /dev/null > $R/run.log 2>&1 ); local rc=$?; kill $mp 2>/dev/null
  echo "ARM $n rc=$rc wall=$((SECONDS-t0))s FATAL=$(grep -c FATAL $R/run.log) NaN=$(tail -n 60 $R/run.log | grep -ciw nan)" \
       "DT_RETRY=$(grep -c '^DT_RETRY cycle=' $R/run.log) peakGiB=$(awk -F, '{v=$2+0; if (v>m) m=v} END {printf "%.1f", m/1024}' $R/mem.txt)" \
       "| $(grep -o 'Picard iterations mean=[^ ]* max=[^ ]* NON-CONVERGED=[^ ]*' $R/run.log | tail -1)" \
       "| $(grep -o 'lagged solves=.*' $R/run.log | tail -1) | $(grep -o 'e_floor cell-solves=.*' $R/run.log | tail -1)" \
       "| $(grep -o 'dt_retry on.*' $R/run.log | head -1) | $(grep 'cycle=' $R/run.log | tail -n 1 | cut -c1-70)"; }
run P  $BP -r $RST $K time/nlim=$N
run B  $BQ -r $RST $K $DR time/nlim=$N
run C  $BQ -r $RST $K $DR $ROL time/nlim=$N
run L1 $BQ -r $RST $K $DR $ROL time/nlim=$NH
run L2 $BQ -r $(lastrst $S/L1) $DR $ROL time/nlim=$N
run F  $BQ -r $RST $K $DR $ROL time/nlim=$N time/dt_retry_force_cycle=$FC time/dt_retry_force_fac=1.0
ok=1
cmpp() { local a=$1 b=$2 tag=$3 r; r=$(python3 $X/cmprst.py $(lastrst $S/$a) $(lastrst $S/$b)); echo "$tag rst: $r"; [[ $r == BITWISE* ]] || ok=0
  local h n; for h in $S/$a/*.hst; do f=$(basename $h); n=$(grep -vc '^#' $S/$b/$f)
    if cmp -s <(grep -v '^#' $h | tail -n $n) <(grep -v '^#' $S/$b/$f); then echo "$tag $f BITWISE (last $n rows)"; else echo "$tag $f DIFF"; ok=0; fi; done
  r=$(BINCONV_DIR=$SCR python3 $SCR/cmplast.py $S/$a $S/$b | tail -1); echo "$tag bins: $r"; [[ $r == *BITWISE* ]] || ok=0; }
cmpp P B  "GATE-BASE    P=B "
cmpp C L2 "GATE-RESTART C=L2"
cmpp C F  "GATE-RETRY   C=F "
[ "$(grep -c '^DT_RETRY cycle=' $S/F/run.log)" = 1 ] || { echo "F: expected exactly 1 DT_RETRY line"; ok=0; }
for a in P B C L1 L2 F; do grep -q FATAL $S/$a/run.log && ok=0; tail -n 60 $S/$a/run.log | grep -qiw nan && ok=0; done
for a in C L2; do o=$(parend $(lastrst $S/$a)); echo "HEADER $a <par_end> offset $o (limit 45056)"; [ "$o" -lt 44000 ] || ok=0; done
echo "seed <par_end> offset $(parend $RST)"
if [ $ok = 1 ]; then echo "GATE PASS $(date +%F_%T)"; date > $S/VERIFIED; cp $S/VERIFIED $D/gate/VERIFIED; else echo "GATE FAIL (do NOT release the ROL2U chain)"; fi
