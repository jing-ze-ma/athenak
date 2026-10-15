#!/bin/bash
# ridort-1014: one IDORT dir, all arms CONCURRENTLY (same binary, same threads):
#   ref   = deck AS IS (FIN5.1 defaults; fin51-r default-off == FIN5.1 21c4839d)
#   r     = deck + rad_m1/implicit_tol_rel=0.1 rad_m1/implicit_tol_rel_max=1e-4                        (option R)
#   o     = r + rad_m1/implicit_opac_lag=pred rad_m1/implicit_opac_lag_passes=2 rad_m1/implicit_gas_newton_u4=true  (ROL2U)
# usage: one.sh <dir> <threads>   env: BIN (fin51-r plain CPU build), [BIN_O: binary for arm o only], ARMS (default "ref r"; set "ref r o" for ROL2U),
#        DECKS (tree with <dir>/in.athinput, args.txt[, two.txt]), DECKS2 (second deck tree, e.g. the bundle's decks_rk),
#        OUT (results root: OUT/runs/<dir>/<arm>), SHOCK (shockref dir for the cart shock decks), HERE (bundle dir), PY, BINCONV
r=$1; nt=$2; s=$DECKS/$r; [ -f $s/in.athinput ] || s=$DECKS2/$r; d=$OUT/runs/$r; mkdir -p $d
KR="rad_m1/implicit_tol_rel=0.1 rad_m1/implicit_tol_rel_max=1e-4"
KO="$KR rad_m1/implicit_opac_lag=pred rad_m1/implicit_opac_lag_passes=2 rad_m1/implicit_gas_newton_u4=true"
arm() { x=$1; { [ -f $d/$x/rc.txt ] && [ "$(cat $d/$x/rc.txt)" = 0 ]; } && return
  b=$BIN; [ $x = o ] && [ -n "$BIN_O" ] && b=$BIN_O
  rm -rf $d/$x; mkdir -p $d/$x; cd $d/$x
  case $x in r) ${PY:-python3} $HERE/setkeys.py $s/in.athinput in.athinput $KR;; o) ${PY:-python3} $HERE/setkeys.py $s/in.athinput in.athinput $KO;;
    *) cp $s/in.athinput in.athinput;; esac
  sed -i -e "s#/raven/ptmp/jinma/suite_1011/shock/#$SHOCK/#" -e "s#/viper/ptmp2/jinma/fin3_1012/shockref/#$SHOCK/#" in.athinput
  t0=$(date +%s); echo "$(md5sum $b)" > bin.txt
  if [ -f $s/two.txt ]; then read nb cb rb < $s/two.txt
    OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $b -i in.athinput -d $d/$x $(cat $s/args.txt) > logA 2>&1
    f=$(ls rst/*.rst | sort | tail -1); xa="time/nlim=$nb rad_m1/implicit_cfl=$cb"; [ -n "$rb" ] && xa="$xa rad_m1/vet_gd_rotate_dt=$rb"
    OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $b -r $f -d $d/$x $xa > run.log 2>&1; echo $? > rc.txt; cat logA run.log > all.log
  else OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $b -i in.athinput -d $d/$x $(cat $s/args.txt) > run.log 2>&1; echo $? > rc.txt; cp run.log all.log; fi
  echo "$(( $(date +%s) - t0 )) s nt=$nt" > time.txt; rm -rf rst; }
for x in ${ARMS:-ref r}; do arm $x & done; wait
cd $d
pic() { grep -h -o 'Picard iterations mean=[^ ]* max=[^ ]* NON-CONVERGED=[^ ]*' $1/all.log | tail -1 | sed 's/Picard iterations //'; }
for x in ${ARMS:-ref r}; do [ $x = ref ] && continue
  echo "$x: $(${PY:-python3} $HERE/cmpout.py $d/ref $d/$x) | rc=$(cat ref/rc.txt)/$(cat $x/rc.txt) fatal=$(grep -c FATAL ref/all.log)/$(grep -c FATAL $x/all.log) dtretry=$(grep -c '^DT_RETRY cycle' ref/all.log)/$(grep -c '^DT_RETRY cycle' $x/all.log) | REF $(pic ref) | $x $(pic $x) | $(cat ref/time.txt)/$(cut -d' ' -f1 $x/time.txt) | $r" > cmp_$x.txt
done
cat cmp_*.txt | tr '\n' ' '; echo
