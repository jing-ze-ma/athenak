#!/bin/bash
# prodg.sh STAR BIN DECK [extra key=val ...]   (orion FRESHGO 2026-11-20: fresh PRODUCTION, Caltech H200, 1 node = 4 GPUs, NR=4)
#   env (required): HOLDER, NODE, OUT (run dir root; the run lives in $OUT/$STAR), FRESH (gunzipped fresh_prod_1120), END (wall end "YYYY-MM-DD HH:MM" local, e.g. holder end - 20 min)
#   env (optional): MEMWRAP, KEEPRST (rst files kept, default 3; older ones deleted after each new one appears)
# Fresh start if $OUT/$STAR/rst has no .rst, else restart from the newest .rst (same command = resubmit). Wall limit -t = END - now - 10 min (final rst written at the limit).
STAR=$1; BIN=$2; DECK=$3; shift 3; EX="$*"
: "${STAR:?}" "${BIN:?}" "${DECK:?}" "${HOLDER:?}" "${NODE:?}" "${OUT:?}" "${FRESH:?}" "${END:?}"; H=$(dirname "$(readlink -f "$0")"); MW=${MEMWRAP:-$H/../memwrap.sh}
D=$OUT/$STAR; mkdir -p $D; cd $D
S=$(( $(date -d "$END" +%s) - $(date +%s) - 600 )); [ $S -gt 1800 ] || { echo "too little wall time ($S s)"; exit 1; }
WT=$(printf "%d:%02d:%02d" $((S/3600)) $((S%3600/60)) $((S%60)))
R=${RST_K9:-$(ls -t $D/rst/*.rst 2>/dev/null | head -1)}; n=$(ls $D/log.* $D/log_k9.* 2>/dev/null | wc -l)
if [ -z "$R" ]; then sed "s#FRESH_DIR#$FRESH#g" $DECK > $D/deck.athinput; grep -q FRESH_DIR $D/deck.athinput && exit 1
  for f in $(sed "s/#.*//" $D/deck.athinput | grep -oE "$FRESH/[^ ]+" | sort -u); do [ -s $f ] || { echo "MISSING $f"; exit 1; }; done; A="-i $D/deck.athinput"
else A="-r $R"; fi
echo "$(date '+%F %T') chunk $n STAR=$STAR BIN=$BIN md5=$(md5sum < $BIN | cut -c1-32) DECK=$DECK md5=$(md5sum < $DECK | cut -c1-32) start=$A WT=$WT HOLDER=$HOLDER NODE=$NODE EX=$EX" >> $D/chunks.txt
( while sleep 600; do ls -t $D/rst/*.rst 2>/dev/null | tail -n +$(( ${KEEPRST:-3} + 1 )) | xargs -r rm -f; done ) & PR=$!
OMPI_MCA_io=romio321 srun --overlap --jobid=$HOLDER -w $NODE -N 1 -n 4 --mpi=pmix -J p_$STAR \
  env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 $MW $D $BIN $A -d $D -t $WT $EX < /dev/null > $D/log_k9.$n 2>&1
echo "$(date '+%F %T') chunk $n rc=$? last: $(grep -E "cycle=" $D/log_k9.$n | tail -1)" >> $D/chunks.txt; kill $PR
