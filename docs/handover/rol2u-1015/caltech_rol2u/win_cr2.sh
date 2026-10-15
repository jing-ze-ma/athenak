#!/bin/bash -l
# win_cr2.sh (CR2): short comparison WINDOW from the SEED copy, tlim = seed t + WIN (hg 2.5e5 s, e1 2e6 s), production output cadence.
#   ARMS="B RW" (default, one job: BASE then ROL2U-window, sequential) or ARMS=H (optional, own job: BASE cfl/2 = dt-error yardstick).
#   sbatch -J cr2w_$CASE -A carnegie_poc -p gpu --qos=normal -t 06:00:00 --nodes=1 --ntasks=$NGPU -c 8 --gres=gpu:nvidia_h200:$NGPU \
#          --mem=240G --exclude=hpc-sm-01-09 --export=ALL,CASE=hg|e1,BQ_MD5=<md5>,CR2_X=<bundle>[,ARMS=H] win_cr2.sh
#   per-arm athena -t (wall cap; the window then ends early and cr2ana.py uses the common part): TB (B, H) / TR (RW), see below.
# B  = 1bd0c2ce production keys + dt_retry (bitwise to the production binary per GATE-BASE);  RW = B + ROL2U (counters for the window:
# the ROL2U chain prints BE steps / floor clips / lagged solves only at the end of a 7-d link);  H = B + time/cfl_number = deck/2.
set -u
source ${CR2_X:?export CR2_X=<this bundle dir>}/cr2_common.sh
chkmd5 $BQ $BQ_MD5
[ -f $D/gate/VERIFIED ] || { echo "gate not passed ($D/gate/VERIFIED missing): exit"; exit 1; }
RST=$(seedrst); read -r _ T0 NC0 _ <<< "$(seedinfo)"
TL=$(python3 -c "print('%.17g' % ($T0 + $WIN))")
CFL=$(grep -a -m1 -E '^ *cfl_number *=' $RST | sed 's/#.*//' | awk -F= '{print $2+0}'); CFL2=$(python3 -c "print($CFL/2)")
case $CASE in hg) TB=${TB:-02:45:00}; TR=${TR:-02:00:00}; TH=${TH:-05:30:00};; e1) TB=${TB:-03:30:00}; TR=${TR:-02:00:00}; TH=${TH:-07:00:00};; esac
K="$(seedkeys) $PKEYS time/tlim=$TL"
echo "CASE $CASE seed t $T0 cycle $NC0 -> tlim $TL; cfl $CFL (H: $CFL2); ARMS ${ARMS:-B RW}; node $(hostname)"
run() { local n=$1 tl=$2; shift 2; local R=$D/win/$n; [ -e $R/run.log ] && { echo "$R exists: skip"; return; }; mkdir -p $R
  echo "$BQ $*" > $R/cmd.txt; local t0=$SECONDS
  ( cd $R; $LAUNCH $BQ -r $RST "$@" -d $R -t $tl < /dev/null > $R/run.log 2>&1 ); local rc=$?
  echo "ARM $n rc=$rc wall=$((SECONDS-t0))s FATAL=$(grep -c FATAL $R/run.log) DT_RETRY=$(grep -c '^DT_RETRY cycle=' $R/run.log)" \
       "| $(grep -o 'Picard iterations mean=[^ ]* max=[^ ]* NON-CONVERGED=[^ ]*' $R/run.log | tail -1) | $(grep 'cycle=' $R/run.log | tail -n 1 | cut -c1-70)"; }
for a in ${ARMS:-B RW}; do case $a in
  B)  run B  $TB $K $DR;;
  RW) run RW $TR $K $DR $ROL;;
  H)  run H  $TH $K $DR time/cfl_number=$CFL2;;
esac; done
