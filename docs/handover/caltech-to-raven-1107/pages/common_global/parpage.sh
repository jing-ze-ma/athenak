#!/bin/bash
# parpage.sh rsgE|hggl [W] (USER 10-22 fast pages): reduce all un-reduced dumps IN PARALLEL (W workers, ONLY=<idx> per worker),
# then update.sh (finds nothing left to reduce: figs + movies + page), then stage to ~/ATHENAK/pages/<name>_glob. Run as a job.
J=/resnick/groups/carnegie_poc/jingze; N=$1; W=${2:-12}
case $N in rsgE) PG=$J/rsgE_page;; hggl) PG=$J/hgglob_page;; *) exit 1;; esac
eval "$(grep -E '^(export|HERE=|R=)' $PG/update.sh)"; export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
cd $PG/out || exit 1
IDX=$(for d in ${BINDIR//:/ }; do ls $d/*.hydro_w.*.bin 2>/dev/null; done | sed -E 's/.*\.([0-9]+)\.bin$/\1/' | sort -u | while read i; do [ -f data/red_$i.npz ] || echo $i; done)
echo "$(date +%T) $N: $(echo $IDX | wc -w) dumps to reduce with $W workers"
echo $IDX | tr ' ' '\n' | grep . | xargs -P $W -I{} env ONLY={} $PY $SRC/reduce.py $RUN data
echo "$(date +%T) reduce done"
DST=/resnick/home/jingze/ATHENAK/pages/${N}_glob
bash $PG/update.sh > /dev/null && rsync -a --delete $PG/out/publish/ $DST/ && echo "$(date '+%m-%d %H:%M') parpage $N" > $DST/STAGED
echo "$(date +%T) staged rc=$?"
