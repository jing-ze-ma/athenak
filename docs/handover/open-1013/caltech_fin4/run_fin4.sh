#!/bin/bash
# fin4-1014 (Raven worker MG): the merged method candidate FIN4 = branch fin4-1014 on the IDORT suite (359 CVALL dirs).
# FIN4 DEFAULTS (fresh run): time2_enth_be = true, implicit_cs2_vimp_npass = 3, implicit_cs2_vimp = true, vet_path = none.
# >>> Run every deck AS IS: do NOT add ANY key (not time2_enth_be, not implicit_cs2_vimp_npass, not implicit_cs2_vimp,
# >>> not vet_path). The test is that FIN4's DEFAULTS reproduce OT's ENTHBE runs.
# bash run_fin4.sh <FIN4 = plain CPU fin4-1014> <OLD = athena_cpu_omp_281dd686 (md5 9ae9f09a)> [parallel=8] [part=A|B|AB]
#  A (list_a.txt, the 214 ENTHBE dirs) in $EB/<dir> (the ENTHBE tree, which has <dir>/enthbe with bin/tab):
#    FIN4 to the end in <dir>/fin4 with the SAME thread count as enthbe/nt.txt (else the list's); two.txt = two-stage
#    (stage A fresh, then -r with time/nlim, implicit_cfl, vet_gd_rotate_dt, as ENTHBE); then
#    cmpout.py enthbe fin4 -> fin4/cmp.txt (BITWISE = every bin after <par_end>, tab and hst data byte-equal).
#  B (list_b.txt, the other 145 dirs) in $MB/<dir> (the mbband tree): BIT, OLD and FIN4 10 cycles each, same threads, final
#    rst -> cmprst_ra.py (the mbband rule: BITWISE / RAHIST = same state; DIFF) -> <dir>/bit4.txt. No FULL runs.
HERE=$(cd $(dirname $0) && pwd); BN=$(readlink -f ${1:?}); BO=$(readlink -f ${2:?}); P=${3:-8}; PART=${4:-AB}; NBIT=10
EB=${EB:-/resnick/groups/carnegie_poc/jingze/enthbe_1013/caltech_enthbe}; MB=${MB:-/resnick/groups/carnegie_poc/jingze/mbband_1013/caltech_mbband}
CRA=${CRA:-$MB/cmprst_ra.py}; [ -f $CRA ] || CRA=$HERE/../../caltech_mbband/cmprst_ra.py
shk() { sed -i -e "s#/raven/ptmp/jinma/suite_1011/shock/#$1/shockref/#" -e "s#/viper/ptmp2/jinma/fin3_1012/shockref/#$1/shockref/#" $2; }
oneA() { r=$1; nt=$2; d=$EB/runs/$r; cd $d || { echo "A NODIR $r"; return; }
  [ -f enthbe/nt.txt ] && nt=$(cat enthbe/nt.txt)
  if ! { [ -f fin4/rc.txt ] && [ "$(cat fin4/rc.txt)" = 0 ]; }; then
    rm -rf fin4; mkdir -p fin4; cd fin4; cp ../in.athinput in.athinput; shk $EB in.athinput; t0=$(date +%s)
    if [ -f ../two.txt ]; then
      read nb cb rb < ../two.txt
      OMP_NUM_THREADS=$nt $BN -i in.athinput -d $d/fin4 $(cat ../args.txt) > logA 2>&1; echo $? > rcA.txt
      f=$(ls rst/*.rst | sort | tail -1); xa="time/nlim=$nb rad_m1/implicit_cfl=$cb"; [ -n "$rb" ] && xa="$xa rad_m1/vet_gd_rotate_dt=$rb"
      OMP_NUM_THREADS=$nt $BN -r $f -d $d/fin4 $xa > run.log 2>&1; echo $? > rc.txt
    else
      OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $BN -i in.athinput -d $d/fin4 $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt
    fi
    echo "$(( $(date +%s) - t0 )) s" > time.txt; echo $nt > nt.txt; rm -rf rst; cd $d
  fi
  echo "$(python3 $HERE/cmpout.py $d/enthbe $d/fin4) $r" > fin4/cmp.txt
  echo "A rc=$(cat fin4/rc.txt) $(cat fin4/time.txt) nt=$nt keys: $(grep -m1 -o 'time2_enth_be=[a-z]* implicit_cs2_vimp_npass=[0-9]*' fin4/run.log) | $(cat fin4/cmp.txt)"; }
oneB() { r=$1; nt=$2; d=$MB/runs/$r; cd $d || { echo "B NODIR $r"; return; }
  if [ ! -f bit4.txt ]; then
    rm -rf bit4O bit4F; mkdir -p bit4O bit4F
    for x in O F; do python3 $HERE/setkeys.py in.athinput bit4$x/in.athinput time/nlim=$NBIT
      grep -qE '^ *file_type *= *rst' bit4$x/in.athinput || printf '\n<output99>\nfile_type = rst\ndt = 1.0e30\n' >> bit4$x/in.athinput
      shk $MB bit4$x/in.athinput; done
    for x in O F; do b=$BN; [ $x = O ] && b=$BO
      ( cd bit4$x; env OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $b -i in.athinput -d $d/bit4$x $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt ) &
    done; wait
    fo=$(ls bit4O/rst/*.rst 2>/dev/null | sort | tail -1); ff=$(ls bit4F/rst/*.rst 2>/dev/null | sort | tail -1)
    if [ -n "$fo" ] && [ -n "$ff" ]; then s=$(python3 $CRA $fo $ff | awk '{print $1}'); else s="FAIL(rc=$(cat bit4O/rc.txt)/$(cat bit4F/rc.txt))"; fi
    echo "$s $r keys: $(grep -m1 -o 'time2_enth_be=[a-z]* implicit_cs2_vimp_npass=[0-9]*' bit4F/run.log)" > bit4.txt; rm -rf bit4O/rst bit4F/rst
  fi
  echo "B $(cat bit4.txt)"; }
export -f oneA oneB shk; export HERE BN BO NBIT EB MB CRA
md5sum $BN $BO
[[ $PART == *B* ]] && { grep -v "^#" $HERE/list_b.txt | xargs -P $P -L 1 bash -c 'oneB "$@"' _; }
[[ $PART == *A* ]] && { grep -v "^#" $HERE/${LISTA:-list_a.txt} | xargs -P $P -L 1 bash -c 'oneA "$@"' _; }
echo "done $(date)"
cat $EB/runs/*/*/*/fin4/cmp.txt $EB/runs/*/*/fin4/cmp.txt 2>/dev/null | sort -u > $HERE/fin4_cmpA.txt; awk '{print $1}' $HERE/fin4_cmpA.txt | sort | uniq -c
cat $MB/runs/*/*/*/bit4.txt $MB/runs/*/*/bit4.txt 2>/dev/null | sort -u > $HERE/fin4_bitB.txt; awk '{print $1}' $HERE/fin4_bitB.txt | sort | uniq -c
