#!/bin/bash
# WV wavefix-1013: the FULL IDORT suite (no bubbles / 2-source) on ONE binary NEW = wavefix-1013 @ 3972a05d
# (fixes 1 implicit_bc_advect_wall + 3 vet_gd_band_dir, both default ON; key 2 implicit_hr_q_face default OFF).
# Arms: (a) fixes 1+3 (deck as is); (b) fixes 1+3 + rad_m1/implicit_hr_q_face=true; (z) keys off (must = OLD).
# bash run_wavefix.sh <OLD = plain CPU 9972a56d> <NEW = plain CPU wavefix-1013 @ 3972a05d> [parallel=16] [NBIT=10]
# Per run dir (list.txt "<dir> <threads>"; in.athinput + args.txt [+ two.txt "nlimB cflB rotB" = two-stage]):
#  1. BIT: OLD, NEW(a), NEW(b), NEW(z), NBIT cycles each, same threads; cmprst.py -> bit.txt
#     "<OLD vs a> <OLD vs b> <OLD vs z> <dir> <s OLD> <s a> <s b> <s z>". (z) MUST be BITWISE everywhere.
#  2. FULL (a) in <dir>/full_a for dirs with (a) DIFF; FULL (b) in <dir>/full_b for dirs with (b) DIFF.
HERE=$(cd $(dirname $0) && pwd); BO=$(readlink -f ${1:?}); BN=$(readlink -f ${2:?}); P=${3:-16}; NBIT=${4:-10}
one() { r=$1; nt=$2; d=$HERE/runs/$r; cd $d || return
  if [ ! -f bit.txt ]; then
    rm -rf bitO bitA bitB bitZ; mkdir -p bitO bitA bitB bitZ
    python3 $HERE/setkeys.py in.athinput bitO/in.athinput time/nlim=$NBIT
    python3 $HERE/setkeys.py in.athinput bitA/in.athinput time/nlim=$NBIT
    python3 $HERE/setkeys.py in.athinput bitB/in.athinput time/nlim=$NBIT rad_m1/implicit_hr_q_face=true
    python3 $HERE/setkeys.py in.athinput bitZ/in.athinput time/nlim=$NBIT rad_m1/implicit_bc_advect_wall=false rad_m1/vet_gd_band_dir=false
    for x in O A B Z; do
      grep -qE '^ *file_type *= *rst' bit$x/in.athinput || printf '\n<output99>\nfile_type = rst\ndt = 1.0e30\n' >> bit$x/in.athinput
      sed -i "s#/raven/ptmp/jinma/suite_1011/shock/#$HERE/shockref/#" bit$x/in.athinput
    done
    for x in O A B Z; do b=$BN; [ $x = O ] && b=$BO
      ( cd bit$x; /usr/bin/time -f "%e" -o sec.txt env OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $b -i in.athinput -d $d/bit$x $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt ) &
    done; wait
    fo=$(ls bitO/rst/*.rst 2>/dev/null | sort | tail -1); fa=$(ls bitA/rst/*.rst 2>/dev/null | sort | tail -1); fz=$(ls bitZ/rst/*.rst 2>/dev/null | sort | tail -1); fb=$(ls bitB/rst/*.rst 2>/dev/null | sort | tail -1)
    if [ -n "$fo" ] && [ -n "$fa" ] && [ -n "$fz" ] && [ -n "$fb" ]; then
      sa=$(python3 $HERE/cmprst.py $fo $fa | awk '{print $1}'); sb=$(python3 $HERE/cmprst.py $fo $fb | awk '{print $1}'); sz=$(python3 $HERE/cmprst.py $fo $fz | awk '{print $1}')
    else sa="FAIL(rc=$(cat bitO/rc.txt)/$(cat bitA/rc.txt)/$(cat bitB/rc.txt)/$(cat bitZ/rc.txt))"; sb=$sa; sz=$sa; fi
    echo "$sa $sb $sz $r $(cat bitO/sec.txt 2>/dev/null | tail -1) $(cat bitA/sec.txt 2>/dev/null | tail -1) $(cat bitB/sec.txt 2>/dev/null | tail -1) $(cat bitZ/sec.txt 2>/dev/null | tail -1)" > bit.txt; rm -rf bitO/rst bitA/rst bitB/rst bitZ/rst
  fi
  echo "$(cat bit.txt)"
  for arm in a b; do col=1; ex=""; [ $arm = b ] && { col=2; ex="rad_m1/implicit_hr_q_face=true"; }
  cd $d
  if [ -z "$NOFULL" ] && { [ "$FULLALL" = 1 ] || [ "$(awk -v c=$col '{print $c}' bit.txt)" != BITWISE ]; } && ! { [ -f full_$arm/rc.txt ] && [ "$(cat full_$arm/rc.txt)" = 0 ]; }; then
    rm -rf full_$arm; mkdir -p full_$arm; cd full_$arm
    python3 $HERE/setkeys.py ../in.athinput in.athinput $ex
    sed -i "s#/raven/ptmp/jinma/suite_1011/shock/#$HERE/shockref/#" in.athinput
    if [ -f ../two.txt ]; then
      read nb cb rb < ../two.txt
      OMP_NUM_THREADS=$nt $BN -i in.athinput -d $d/full_$arm $(cat ../args.txt) > logA 2>&1
      f=$(ls rst/*.rst | sort | tail -1)
      xa="time/nlim=$nb rad_m1/implicit_cfl=$cb"; [ -n "$rb" ] && xa="$xa rad_m1/vet_gd_rotate_dt=$rb"
      /usr/bin/time -f "%e s" -o time.txt env OMP_NUM_THREADS=$nt $BN -r $f -d $d/full_$arm $xa > run.log 2>&1; echo $? > rc.txt
    else
      /usr/bin/time -f "%e s" -o time.txt env OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $BN -i in.athinput -d $d/full_$arm $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt
    fi
    rm -rf rst; echo "FULL_$arm $r rc=$(cat rc.txt) $(cat time.txt) $(grep -o 'Picard iterations mean=[0-9.e+-]*' run.log | tail -1) $(grep -o 'reach theta [0-9]* phi [0-9]*' run.log | head -1) $(grep -o 'lateral reads clamped[^<]*' run.log | tail -1)"
  fi
  done; }
export -f one; export HERE BO BN NBIT NOFULL FULLALL
md5sum $BO $BN; xargs -P $P -L 1 bash -c 'one "$@"' _ < $HERE/list.txt; echo "done $(date)"
cat $HERE/runs/*/*/*/bit.txt $HERE/runs/ss/*/bit.txt 2>/dev/null | sort -k3 > $HERE/bit_all.txt
echo "arm a (fixes 1+3) vs OLD:"; awk '{print $1}' $HERE/bit_all.txt | sort | uniq -c; echo "arm b (+ key 2) vs OLD:"; awk '{print $2}' $HERE/bit_all.txt | sort | uniq -c; echo "arm z (keys off) vs OLD (must be all BITWISE):"; awk '{print $3}' $HERE/bit_all.txt | sort | uniq -c
