#!/bin/bash
# movies_par.sh rsgE|hggl [NSHARD] (10-23): FULL parallel re-render of the 2-D movie frames (dpi change etc.): limits pass,
# clear frames, NSHARD render processes, encode, stage movies + fig_eq_latest to the page dir. Run in an expansion job.
J=/resnick/groups/carnegie_poc/jingze; N=$1; NS=${2:-16}
case $N in rsgE) PG=$J/rsgE_page;; hggl) PG=$J/hgglob_page;; esac
eval "$(grep -E '^(export|HERE=|R=)' $PG/update.sh)"; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
cd $PG/out; exec 9>$PG/out/.update.lock; flock 9
LIMITS_ONLY=1 $PY $SRC/mk_movies.py data $PG/out
rm -f frames/{eq,inner,phot,map}/*.png
for k in $(seq 0 $((NS-1))); do SHARD=$k/$NS $PY $SRC/mk_movies.py data $PG/out > frames/shard$k.log 2>&1 & done; wait
$PY $SRC/mk_movies.py data $PG/out && rsync -a --exclude "*_3d*" $PG/out/publish/movies/ /resnick/home/jingze/ATHENAK/pages/${N}_glob/movies/
echo "$N movies done $(date +%T)"
