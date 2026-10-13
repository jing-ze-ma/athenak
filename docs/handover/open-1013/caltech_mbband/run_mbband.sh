#!/bin/bash
# OD mbband-1013: the FULL IDORT suite (no bubbles / 2-source) in two arms on ONE binary NEW = mbband-1013 @ 34f79349:
#   (a) NEW, deck as is (vet_path default none): band fix alone; (b) NEW + rad_m1/vet_path=cells (vet_sc and vet_gd).
# bash run_mbband.sh <OLD = plain CPU 281dd686> <NEW = plain CPU mbband-1013 @ 34f79349> [parallel=16] [NBIT=10]
# Per run dir (list.txt "<dir> <threads>"; in.athinput + args.txt [+ two.txt "nlimB cflB rotB" = two-stage]):
#  1. BIT: OLD, NEW(a), NEW(b) NBIT cycles each, same threads, final rst; cmprst.py -> bit.txt:
#     "<OLD vs a> <OLD vs b> <dir> <sec OLD> <sec a> <sec b>"  (a must be BITWISE everywhere).
#  2. FULL (b) for EVERY dir (cells can switch on after NBIT cycles): NEW + vet_path=cells to the end, in <dir>/full_b.
HERE=$(cd $(dirname $0) && pwd); BO=$(readlink -f ${1:?}); BN=$(readlink -f ${2:?}); P=${3:-16}; NBIT=${4:-10}
one() { r=$1; nt=$2; d=$HERE/runs/$r; cd $d || return
  if [ ! -f bit.txt ]; then
    rm -rf bitO bitA bitB; mkdir -p bitO bitA bitB
    python3 $HERE/setkeys.py in.athinput bitO/in.athinput time/nlim=$NBIT
    python3 $HERE/setkeys.py in.athinput bitA/in.athinput time/nlim=$NBIT
    python3 $HERE/setkeys.py in.athinput bitB/in.athinput time/nlim=$NBIT rad_m1/vet_path=cells
    for x in O A B; do
      grep -qE '^ *file_type *= *rst' bit$x/in.athinput || printf '\n<output99>\nfile_type = rst\ndt = 1.0e30\n' >> bit$x/in.athinput
      sed -i "s#/raven/ptmp/jinma/suite_1011/shock/#$HERE/shockref/#" bit$x/in.athinput
    done
    for x in O A B; do b=$BN; [ $x = O ] && b=$BO
      ( cd bit$x; /usr/bin/time -f "%e" -o sec.txt env OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $b -i in.athinput -d $d/bit$x $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt ) &
    done; wait
    fo=$(ls bitO/rst/*.rst 2>/dev/null | sort | tail -1); fa=$(ls bitA/rst/*.rst 2>/dev/null | sort | tail -1); fb=$(ls bitB/rst/*.rst 2>/dev/null | sort | tail -1)
    if [ -n "$fo" ] && [ -n "$fa" ] && [ -n "$fb" ]; then
      sa=$(python3 $HERE/cmprst.py $fo $fa | awk '{print $1}'); sb=$(python3 $HERE/cmprst.py $fo $fb | awk '{print $1}')
    else sa="FAIL(rc=$(cat bitO/rc.txt)/$(cat bitA/rc.txt)/$(cat bitB/rc.txt))"; sb=$sa; fi
    echo "$sa $sb $r $(cat bitO/sec.txt 2>/dev/null | tail -1) $(cat bitA/sec.txt 2>/dev/null | tail -1) $(cat bitB/sec.txt 2>/dev/null | tail -1)" > bit.txt; rm -rf bitO/rst bitA/rst bitB/rst
  fi
  echo "$(cat bit.txt)"
  if [ -z "$NOFULL" ] && ! { [ -f full_b/rc.txt ] && [ "$(cat full_b/rc.txt)" = 0 ]; }; then
    rm -rf full_b; mkdir -p full_b; cd full_b
    python3 $HERE/setkeys.py ../in.athinput in.athinput rad_m1/vet_path=cells
    sed -i "s#/raven/ptmp/jinma/suite_1011/shock/#$HERE/shockref/#" in.athinput
    if [ -f ../two.txt ]; then
      read nb cb rb < ../two.txt
      OMP_NUM_THREADS=$nt $BN -i in.athinput -d $d/full_b $(cat ../args.txt) > logA 2>&1
      f=$(ls rst/*.rst | sort | tail -1)
      xa="time/nlim=$nb rad_m1/implicit_cfl=$cb"; [ -n "$rb" ] && xa="$xa rad_m1/vet_gd_rotate_dt=$rb"
      /usr/bin/time -f "%e s" -o time.txt env OMP_NUM_THREADS=$nt $BN -r $f -d $d/full_b $xa > run.log 2>&1; echo $? > rc.txt
    else
      /usr/bin/time -f "%e s" -o time.txt env OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $BN -i in.athinput -d $d/full_b $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt
    fi
    rm -rf rst; echo "FULL_b $r rc=$(cat rc.txt) $(cat time.txt) $(grep -o 'subdivided steps (ALL ranks, all sweeps) [0-9.e+]*' run.log | tail -1)"
  fi; }
export -f one; export HERE BO BN NBIT NOFULL
md5sum $BO $BN; xargs -P $P -L 1 bash -c 'one "$@"' _ < $HERE/list.txt; echo "done $(date)"
cat $HERE/runs/*/*/*/bit.txt $HERE/runs/ss/*/bit.txt 2>/dev/null | sort -k3 > $HERE/bit_all.txt
echo "arm a:"; awk '{print $1}' $HERE/bit_all.txt | sort | uniq -c; echo "arm b:"; awk '{print $2}' $HERE/bit_all.txt | sort | uniq -c
