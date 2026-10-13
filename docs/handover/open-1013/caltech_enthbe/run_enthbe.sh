#!/bin/bash
# open-ot-1013 ENTHBE: the candidate time2_enth_be = true + implicit_cs2_vimp_npass = 3 (+ implicit_cs2_vimp = true) on
# every active IDORT run (the 214 dirs that were DIFF in the CVALL BIT step; Raven checks the other 145 bitwise).
# bash run_enthbe.sh <NEW binary = 281dd686 + patch/0001-0006 = open-ot-1013-ht @ c685e6dd> [parallel=8]
# Per dir: NEW, deck + the three keys (named), to the end -> <dir>/enthbe (tab/ bin/ run.log rc.txt time.txt).
# two.txt = two-stage (patm*): stage A, then -r with time/nlim, implicit_cfl, vet_gd_rotate_dt.
HERE=$(cd $(dirname $0) && pwd); BN=$(readlink -f ${1:?}); P=${2:-8}
one() { r=$1; nt=$2; d=$HERE/runs/$r; cd $d || return
  [ -f enthbe/rc.txt ] && [ "$(cat enthbe/rc.txt)" = 0 ] && return
  rm -rf enthbe; mkdir -p enthbe; cd enthbe
  python3 $HERE/setkeys.py ../in.athinput in.athinput rad_m1/time2_enth_be=true rad_m1/implicit_cs2_vimp_npass=3 rad_m1/implicit_cs2_vimp=true
  sed -i -e "s#/raven/ptmp/jinma/suite_1011/shock/#$HERE/shockref/#" -e "s#/viper/ptmp2/jinma/fin3_1012/shockref/#$HERE/shockref/#" in.athinput
  t0=$(date +%s)
  if [ -f ../two.txt ]; then
    read nb cb rb < ../two.txt
    OMP_NUM_THREADS=$nt $BN -i in.athinput -d $d/enthbe $(cat ../args.txt) > logA 2>&1
    f=$(ls rst/*.rst | sort | tail -1); xa="time/nlim=$nb rad_m1/implicit_cfl=$cb"; [ -n "$rb" ] && xa="$xa rad_m1/vet_gd_rotate_dt=$rb"
    OMP_NUM_THREADS=$nt $BN -r $f -d $d/enthbe $xa > run.log 2>&1; echo $? > rc.txt
  else
    OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $BN -i in.athinput -d $d/enthbe $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt
  fi
  echo "$(( $(date +%s) - t0 )) s" > time.txt; rm -rf rst; echo "FULL $r rc=$(cat rc.txt) $(cat time.txt)"; }
export -f one; export HERE BN
md5sum $BN; grep -v "^#" $HERE/${LIST:-list.txt} | xargs -P $P -L 1 bash -c 'one "$@"' _; echo "done $(date)"
