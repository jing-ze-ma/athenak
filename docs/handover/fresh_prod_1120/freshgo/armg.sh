#!/bin/bash
# armg.sh TAG BIN DECK NC [EXTRA key=val ...]   (orion FRESHGO 2026-11-20: fresh-start smoke / block-shape timing, Caltech H200, 1 node)
#   env (required): HOLDER, NODE, OUT, FRESH (abs path of a GUNZIPPED copy of bsg-files docs/handover/fresh_prod_1120; replaces FRESH_DIR in DECK)
#   env (optional): NR (default 4), TL (srun minutes, default 60), RT (1 = also do the restart arm b, default 1), KEEPBIN, MEMWRAP (default $H/../memwrap.sh)
# a  = fresh start -> NC cycles straight;  b = fresh start -> NC/2 (rst) then restart -> NC;  resf.py compares hst rows NC/2+1..NC + cycle-NC bins (bitwise expected)
TAG=$1; BIN=$2; DECK=$3; NC=$4; shift 4; EX="$*"
: "${TAG:?}" "${BIN:?}" "${DECK:?}" "${NC:?}" "${HOLDER:?}" "${NODE:?}" "${OUT:?}" "${FRESH:?}"
case "$TAG" in */*|.|..) echo "bad TAG $TAG"; exit 1;; esac; [ $(echo "$OUT" | tr -cd / | wc -c) -ge 2 ] || { echo "OUT $OUT too shallow"; exit 1; }
NR=${NR:-4}; H=$(dirname "$(readlink -f "$0")"); MW=${MEMWRAP:-$H/../memwrap.sh}
D=$OUT/$TAG; rm -rf $D; mkdir -p $D; cd $D
sed "s#FRESH_DIR#$FRESH#g" $DECK > $D/deck.athinput; grep -q FRESH_DIR $D/deck.athinput && { echo "FRESH_DIR left" > res.txt; exit 1; }
for f in $(grep -oE "$FRESH/[^ #]+" $D/deck.athinput | sort -u); do [ -s $f ] || { echo "MISSING $f" >> res.txt; exit 1; }; done
echo "TAG=$TAG BIN=$BIN md5=$(md5sum < $BIN | cut -c1-32) DECK=$DECK md5=$(md5sum < $DECK | cut -c1-32) NC=$NC HOLDER=$HOLDER NODE=$NODE NR=$NR EX=$EX" > info.txt
run() { local d=$1 n=$2; shift 2; mkdir -p $D/$d; cd $D/$d; t0=$(date +%s)
  OMPI_MCA_io=romio321 srun --overlap --jobid=$HOLDER -w $NODE -N 1 -n $NR --mpi=pmix -t ${TL:-60} -J g_$TAG \
    env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 $MW $D $BIN "$@" time/nlim=$n $EX -d $D/$d \
    < /dev/null > $([ -f run.log ] && echo run2.log || echo run.log) 2>&1
  echo "$d rc=$? wall=$(( $(date +%s)-t0 ))" >> $D/$d/rc.txt; cd $D; }
run a $NC -i $D/deck.athinput
if [ "${RT:-1}" = 1 ]; then H2=$((NC/2)); run b $H2 -i $D/deck.athinput; f=$(ls -t $D/b/rst/*.rst | head -1); echo "restart from $f" >> $D/b/rc.txt; run b $NC -r $f; fi
find $D -name "*.rst" -exec ls -l {} + > $D/rstsizes.txt 2>/dev/null; find $D -name "*.rst" -delete
python3 $H/resf.py $D $((NC/2)) > $D/res.txt 2>&1; cat $D/*/rc.txt >> $D/res.txt
[ "${KEEPBIN:-0}" = 1 ] || find $D -name "*.bin" -delete
echo DONE >> $D/res.txt
