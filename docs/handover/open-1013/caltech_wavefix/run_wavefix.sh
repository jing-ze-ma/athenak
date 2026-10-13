#!/bin/bash
# WV wavefix-1013: the FULL IDORT suite (no bubbles / 2-source) on ONE binary NEW = wavefix-1013 @ 41c8b210
# (fixes 1 implicit_bc_advect_wall + 3 vet_gd_band_dir, both default ON; key 2 implicit_hr_q_face stays OFF).
# bash run_wavefix.sh <OLD = plain CPU 9972a56d> <NEW = plain CPU wavefix-1013 @ 41c8b210> [parallel=16] [NBIT=10]
# Per run dir (list.txt "<dir> <threads>"; in.athinput + args.txt [+ two.txt "nlimB cflB rotB" = two-stage]):
#  1. BIT: OLD, NEW(a) (deck as is: fixes 1+3 on), NEW(z) (+ implicit_bc_advect_wall=false vet_gd_band_dir=false),
#     NBIT cycles each, same threads; cmprst.py -> bit.txt "<OLD vs a> <OLD vs z> <dir> <s OLD> <s a> <s z>".
#     (z) MUST be BITWISE everywhere (keys off = 9972a56d). (a) DIFF where a fix acts.
#  2. FULL (a): NEW to the end in <dir>/full_a, ONLY for dirs with (a) DIFF (or when FULLALL=1).
HERE=$(cd $(dirname $0) && pwd); BO=$(readlink -f ${1:?}); BN=$(readlink -f ${2:?}); P=${3:-16}; NBIT=${4:-10}
one() { r=$1; nt=$2; d=$HERE/runs/$r; cd $d || return
  if [ ! -f bit.txt ]; then
    rm -rf bitO bitA bitZ; mkdir -p bitO bitA bitZ
    python3 $HERE/setkeys.py in.athinput bitO/in.athinput time/nlim=$NBIT
    python3 $HERE/setkeys.py in.athinput bitA/in.athinput time/nlim=$NBIT
    python3 $HERE/setkeys.py in.athinput bitZ/in.athinput time/nlim=$NBIT rad_m1/implicit_bc_advect_wall=false rad_m1/vet_gd_band_dir=false
    for x in O A Z; do
      grep -qE '^ *file_type *= *rst' bit$x/in.athinput || printf '\n<output99>\nfile_type = rst\ndt = 1.0e30\n' >> bit$x/in.athinput
      sed -i "s#/raven/ptmp/jinma/suite_1011/shock/#$HERE/shockref/#" bit$x/in.athinput
    done
    for x in O A Z; do b=$BN; [ $x = O ] && b=$BO
      ( cd bit$x; /usr/bin/time -f "%e" -o sec.txt env OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $b -i in.athinput -d $d/bit$x $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt ) &
    done; wait
    fo=$(ls bitO/rst/*.rst 2>/dev/null | sort | tail -1); fa=$(ls bitA/rst/*.rst 2>/dev/null | sort | tail -1); fz=$(ls bitZ/rst/*.rst 2>/dev/null | sort | tail -1)
    if [ -n "$fo" ] && [ -n "$fa" ] && [ -n "$fz" ]; then
      sa=$(python3 $HERE/cmprst.py $fo $fa | awk '{print $1}'); sz=$(python3 $HERE/cmprst.py $fo $fz | awk '{print $1}')
    else sa="FAIL(rc=$(cat bitO/rc.txt)/$(cat bitA/rc.txt)/$(cat bitZ/rc.txt))"; sz=$sa; fi
    echo "$sa $sz $r $(cat bitO/sec.txt 2>/dev/null | tail -1) $(cat bitA/sec.txt 2>/dev/null | tail -1) $(cat bitZ/sec.txt 2>/dev/null | tail -1)" > bit.txt; rm -rf bitO/rst bitA/rst bitZ/rst
  fi
  echo "$(cat bit.txt)"
  if [ -z "$NOFULL" ] && { [ "$FULLALL" = 1 ] || [ "$(awk '{print $1}' bit.txt)" != BITWISE ]; } && ! { [ -f full_a/rc.txt ] && [ "$(cat full_a/rc.txt)" = 0 ]; }; then
    rm -rf full_a; mkdir -p full_a; cd full_a
    cp ../in.athinput in.athinput
    sed -i "s#/raven/ptmp/jinma/suite_1011/shock/#$HERE/shockref/#" in.athinput
    if [ -f ../two.txt ]; then
      read nb cb rb < ../two.txt
      OMP_NUM_THREADS=$nt $BN -i in.athinput -d $d/full_a $(cat ../args.txt) > logA 2>&1
      f=$(ls rst/*.rst | sort | tail -1)
      xa="time/nlim=$nb rad_m1/implicit_cfl=$cb"; [ -n "$rb" ] && xa="$xa rad_m1/vet_gd_rotate_dt=$rb"
      /usr/bin/time -f "%e s" -o time.txt env OMP_NUM_THREADS=$nt $BN -r $f -d $d/full_a $xa > run.log 2>&1; echo $? > rc.txt
    else
      /usr/bin/time -f "%e s" -o time.txt env OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $BN -i in.athinput -d $d/full_a $(cat ../args.txt) > run.log 2>&1; echo $? > rc.txt
    fi
    rm -rf rst; echo "FULL_a $r rc=$(cat rc.txt) $(cat time.txt) $(grep -o 'Picard iterations mean=[0-9.e+-]*' run.log | tail -1) $(grep -o 'reach theta [0-9]* phi [0-9]*' run.log | head -1) $(grep -o 'lateral reads clamped[^<]*' run.log | tail -1)"
  fi; }
export -f one; export HERE BO BN NBIT NOFULL FULLALL
md5sum $BO $BN; xargs -P $P -L 1 bash -c 'one "$@"' _ < $HERE/list.txt; echo "done $(date)"
cat $HERE/runs/*/*/*/bit.txt $HERE/runs/ss/*/bit.txt 2>/dev/null | sort -k3 > $HERE/bit_all.txt
echo "arm a (fixes on) vs OLD:"; awk '{print $1}' $HERE/bit_all.txt | sort | uniq -c; echo "arm z (keys off) vs OLD (must be all BITWISE):"; awk '{print $2}' $HERE/bit_all.txt | sort | uniq -c
