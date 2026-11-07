#!/bin/bash
#SBATCH -A carnegie_poc -p expansion -N 1 -c 32 --mem=200G -t 01:30:00 -J pgcatch -o /resnick/groups/carnegie_poc/jingze/catchup_job.%j.out
# catchup_job.sh rsgE|hggl (10-24, USER "make it fast"): quickpage catch-up on an expansion node: reduce new dumps 24-wide,
# 2-D frames in NSH (default 24; 10-30: rsgE OOM at 24 x 2128 dumps -> NSH=8) render shards (limits pass first, so a limits redo is parallel too), encode, stage; then 3-D (r3d_par.sh)
# pagefast-1102 (USER "make the refresh pipelines very fast"; old version catchup_job.sh.pre_fast_1102), same outputs:
#  - reduce: new dumps from 2 listings (no per-dump stat); RED_THREADS EOS threads per dump when there are few new dumps
#    (32 cores = processes x threads); m1_vet header times cached (reduce.py)
#  - 2-D movies: ONE parallel limits pass (POOL_PROCS 24) writes limits + redo flags + the missing-frame task list; the shards and
#    the final encode pass read them (LIMITS_FROM=1) instead of re-pooling all ~3900 dumps each (that was ~90 % of the time);
#    shards no longer hold the pools, so NSH 24 fits in memory again; a limits redo now re-encodes that movie in full
#  - per-stage timing lines in the out file
J=/resnick/groups/carnegie_poc/jingze; N=$1
case $N in rsgE) PG=$J/rsgE_page;; hggl) PG=$J/hgglob_page;; esac
eval "$(grep -E '^(export|HERE=|R=)' $PG/update.sh)"; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
DST=/resnick/home/jingze/ATHENAK/pages/${N}_glob
exec 7>$PG/out/.catchup.lock; flock 7; exec 9>$PG/out/.update.lock; flock 9
cd $PG/out
T0=$SECONDS; echo "$N start $(date +%T)"
NEW=$(comm -23 <(ls $BINDIR | sed -nE 's/.*\.hydro_w\.([0-9]+)\.bin$/\1/p' | sort) <(ls data | sed -nE 's/^red_([0-9]+)\.npz$/\1/p' | sort))
K=$(echo -n "$NEW" | grep -c .)
TH=$(( K > 0 ? 32 / K : 1 )); [ $TH -gt 6 ] && TH=6; [ $TH -lt 1 ] && TH=1
NP=$(( 32 / TH )); [ $NP -gt 24 ] && NP=24
echo -n "$NEW" | xargs -r -P $NP -I{} env ONLY={} RED_THREADS=$TH $PY $SRC/reduce.py $RUN data > /dev/null 2>&1
echo "$N reduced $(ls data | grep -c '^red_.*\.npz$') ($K new, $NP x $TH threads, $((SECONDS-T0)) s) $(date +%T)"
T1=$SECONDS
{ LIMITS_ONLY=1 POOL_PROCS=24 $PY $SRC/mk_movies.py data $PG/out && echo "limits pass $((SECONDS-T1)) s" && \
  { NS=${NSH:-24}; for k in $(seq 0 $((NS-1))); do LIMITS_FROM=1 SHARD=$k/$NS $PY $SRC/mk_movies.py data $PG/out > /dev/null 2>&1 & done; wait
    echo "frames $((SECONDS-T1)) s"; LIMITS_FROM=1 $PY $SRC/mk_movies.py data $PG/out; echo "encode $((SECONDS-T1)) s"; }; } > $PG/out/catchup_movies.log 2>&1 \
  && rsync -a $PG/out/publish/movies/ $DST/movies/ --exclude "*_3d*"
echo "$N 2-D movies rc=$? ($((SECONDS-T1)) s: $(grep -E '^(limits pass|frames|encode) ' $PG/out/catchup_movies.log | paste -sd ' ')) $(date +%T)"
flock -u 9
T1=$SECONDS; bash $J/r3d_par.sh $N; echo "$N 3-D rc=$? ($((SECONDS-T1)) s) $(date +%T)"
# 11-01: rebuild the page text from all reduced dumps, then mark READY with the files to publish
T1=$SECONDS; bash $J/quickpage.sh $N nocatch > /dev/null 2>&1; echo "$N page text $((SECONDS-T1)) s"
# movfix-1102: 2-D movies are parts <prefix>_<movie>_pNN.mp4 (<= 14 MB); READY lists only the parts this catch-up wrote
# (frames/changed_2d.txt; closed parts never change) + the 3-D files; the old single-file 2-D movies are gone from staging
rm -f $DST/movies/${N}_eq.mp4 $DST/movies/${N}_inner.mp4 $DST/movies/${N}_phot.mp4 $DST/movies/${N}_map.mp4
( cd $DST && echo "READY $(date '+%m-%d %H:%M') $(grep -oE 't = [0-9.e+]+ s = <b>[0-9.]+ d' index.html | sed 's/<b>//' | head -1)"; cat $PG/out/frames/changed_2d.txt 2>/dev/null; ls movies/*_3d*.mp4 movies/*_poster.png ) > $DST/READY
echo "$N catch-up total $((SECONDS-T0)) s $(date +%T)"
