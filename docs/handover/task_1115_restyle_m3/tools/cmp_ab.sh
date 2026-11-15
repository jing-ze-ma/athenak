#!/bin/bash
# cmp_ab.sh <pkg> <bindir> <name> <prob> <deck> <keys> <links>: run the deck with
# <bindir>/athena_<prob>_A and _B (1 rank) in <pkg>/runs/<name>/{A,B}; compare every
# output file (bincmp.py: data bytes; the embedded parameter text may differ by removed
# keys) and the cycle lines of stdout.  Prints "<name> IDENTICAL|DIFFER|FAIL ...".
PKG=$1; BD=$2; NAME=$3; PROB=$4; DECK=${5//@PKG@/$PKG}; KEYS=${6//@PKG@/$PKG}
LINKS=${7//@PKG@/$PKG}
D=$PKG/runs/$NAME; mkdir -p $D/A $D/B
python3 $PKG/tools/mkdeck.py $DECK $D/deck.athinput $KEYS || exit 1
export OMP_NUM_THREADS=1
for x in A B; do
  cd $D/$x || exit 1
  find . -mindepth 1 -delete
  for l in $LINKS; do ln -s ${l#*=} ${l%%=*}; done
  timeout ${TMO:-7200} $BD/athena_${PROB}_$x -i $D/deck.athinput -d $D/$x > stdout.txt 2>&1
  echo $? > rc.txt
done
cd $D || exit 1
for x in A B; do grep -o 'cycle=.*' $x/stdout.txt | grep -v cpu_second > $x/cycles.txt; done
nf=0; nd=0; diffs=""
for f in $(cd A && find . -type f ! -name stdout.txt ! -name rc.txt ! -name '*.rst' | sort); do
  nf=$((nf+1))
  python3 $PKG/tools/bincmp.py A/$f B/$f 2>/dev/null || { nd=$((nd+1)); diffs="$diffs $f"; }
done
nb=$(cd B && find . -type f ! -name stdout.txt ! -name rc.txt ! -name '*.rst' | wc -l)
ra=$(cat A/rc.txt); rb=$(cat B/rc.txt); ncyc=$(wc -l < A/cycles.txt)
if [ "$ra" != 0 ] || [ "$rb" != 0 ]; then st=FAIL
elif [ $nd = 0 ] && [ $nf = $nb ] && [ $ncyc -gt 0 ]; then st=IDENTICAL
else st=DIFFER; fi
echo "$NAME $st rc=$ra/$rb files=$nf/$nb cycles=$ncyc ndiff=$nd${diffs:+ :$diffs}"
