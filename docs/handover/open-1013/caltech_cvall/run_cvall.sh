#!/bin/bash
# open-ot-1013 CVALL: implicit_cs2_vimp = false across the FULL IDORT suite (no bubbles / 2-source).
# bash run_cvall.sh <OLD binary = plain CPU 281dd686> <NEW binary = 281dd686 + patch/0001> [parallel=16] [NBIT=10]
# Per run dir (list.txt "<dir> <threads>"; in.athinput + args.txt [+ two.txt "nlimB cflB rotB" = two-stage]):
#  1. BIT: OLD (deck as is) and NEW (deck + implicit_cs2_vimp = false), NBIT cycles each, same threads, final rst;
#     cmprst.py compares the data after <par_end> -> bit.txt line "BITWISE|DIFF <dir>".
#  2. FULL: only when DIFF and the dir is not in full_skip.txt (those are in the caltech_cv bundle): NEW to the end
#     (two-stage: stage A, then -r with time/nlim=A+B, implicit_cfl, vet_gd_rotate_dt), outputs in <dir>/full.
HERE=$(cd $(dirname $0) && pwd); BO=$(readlink -f ${1:?}); BN=$(readlink -f ${2:?}); P=${3:-16}; NBIT=${4:-10}
one() { r=$1; nt=$2; d=$HERE/runs/$r; cd $d || return
  if [ ! -f bit.txt ]; then
    rm -rf bitO bitN; mkdir -p bitO bitN
    python3 $HERE/setkeys.py in.athinput bitO/in.athinput time/nlim=$NBIT
    python3 $HERE/setkeys.py in.athinput bitN/in.athinput time/nlim=$NBIT rad_m1/implicit_cs2_vimp=false
    for x in O N; do printf '\n<output99>\nfile_type = rst\ndt = 1.0e30\n' >> bit$x/in.athinput; done
    ( cd bitO; OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $BO -i in.athinput -d $d/bitO $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt ) &
    ( cd bitN; OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $BN -i in.athinput -d $d/bitN $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt ) &
    wait
    fo=$(ls bitO/rst/*.rst 2>/dev/null | sort | tail -1); fn=$(ls bitN/rst/*.rst 2>/dev/null | sort | tail -1)
    if [ -n "$fo" ] && [ -n "$fn" ]; then s=$(python3 $HERE/cmprst.py $fo $fn | awk '{print $1}'); else s="FAIL(rcO=$(cat bitO/rc.txt),rcN=$(cat bitN/rc.txt))"; fi
    echo "$s $r" > bit.txt; rm -rf bitO/rst bitN/rst
  fi
  echo "$(cat bit.txt)"
  if grep -q '^DIFF' bit.txt && ! grep -qx "$r" $HERE/full_skip.txt && ! { [ -f full/rc.txt ] && [ "$(cat full/rc.txt)" = 0 ]; }; then
    rm -rf full; mkdir -p full; cd full
    python3 $HERE/setkeys.py ../in.athinput in.athinput rad_m1/implicit_cs2_vimp=false
    if [ -f ../two.txt ]; then
      read nb cb rb < ../two.txt
      OMP_NUM_THREADS=$nt $BN -i in.athinput -d $d/full $(cat ../args.txt) > logA 2>&1
      f=$(ls rst/*.rst | sort | tail -1)
      xa="time/nlim=$nb rad_m1/implicit_cfl=$cb"; [ -n "$rb" ] && xa="$xa rad_m1/vet_gd_rotate_dt=$rb"
      /usr/bin/time -f "%e s" -o time.txt env OMP_NUM_THREADS=$nt $BN -r $f -d $d/full $xa > run.log 2>&1; echo $? > rc.txt
    else
      /usr/bin/time -f "%e s" -o time.txt env OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $BN -i in.athinput -d $d/full $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt
    fi
    rm -rf rst; echo "FULL $r rc=$(cat rc.txt) $(cat time.txt)"
  fi; }
export -f one; export HERE BO BN NBIT
md5sum $BO $BN; xargs -P $P -L 1 bash -c 'one "$@"' _ < $HERE/list.txt; echo "done $(date)"
cat $HERE/runs/*/*/*/bit.txt $HERE/runs/ss/*/bit.txt 2>/dev/null | sort > $HERE/bit_all.txt; awk '{print $1}' $HERE/bit_all.txt | sort | uniq -c
