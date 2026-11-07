#!/bin/bash
# quickpage.sh rsgE|hggl (USER 10-22 "url updates must be really fast"): login-node page rebuild from hst + logs + the dumps
# ALREADY reduced (no reduce.py, no movies, no 3-D), ~15 s. Then publish index.html (see memory global-pages-sync).
# 10-23 fix: reuses the last fig_eq_latest.jpg (was "missing figure"), and starts a detached CATCH-UP (reduce new dumps 4-wide
# + incremental 2-D frames/movies, flock) so the NEXT quickpage is current. "catchup" arg = the catch-up itself.
J=/resnick/groups/carnegie_poc/jingze; N=$1
case $N in rsgE) PG=$J/rsgE_page;; hggl) PG=$J/hgglob_page;; *) echo "usage: quickpage.sh rsgE|hggl [catchup|status]"; exit 1;; esac
[ "$2" = status ] && { D=/resnick/home/jingze/ATHENAK/pages/${N}_glob; echo "$N: catch-up $(flock -n $PG/out/.catchup.lock true && echo idle || echo BUSY); $(head -1 $D/READY 2>/dev/null || echo "not ready")"; exit 0; }   # 11-01 instant status (before the slow sync)
eval "$(grep -E '^(export|HERE=|R=)' $PG/update.sh)"; export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
DST=/resnick/home/jingze/ATHENAK/pages/${N}_glob; CL=$PG/out/catchup.log
new() { comm -23 <(ls $BINDIR | sed -nE 's/.*\.hydro_w\.([0-9]+)\.bin$/\1/p' | sort) <(ls $PG/out/data | sed -nE 's/^red_([0-9]+)\.npz$/\1/p' | sort); }   # 11-01: 2 listings, no per-file stat (was ~1 min on NFS)
if [ "$2" = status ]; then echo "$N: catch-up $(flock -n $PG/out/.catchup.lock true && echo idle || echo BUSY); $(head -1 $DST/READY 2>/dev/null || echo "not ready")"; exit 0; fi
if [ "$2" = catchup ]; then
  exec 7>$PG/out/.catchup.lock; flock -n 7 || exit 0; exec 9>$PG/out/.update.lock; flock 9   # waits for a running update.sh
  { echo "$(date +%T) catch-up start"; ulimit -v 9000000; cd $PG/out
    new | xargs -r -P 4 -I{} env ONLY={} nice -n 10 $PY $SRC/reduce.py $RUN data > /dev/null 2>&1
    nice -n 10 $PY $SRC/mk_movies.py data $PG/out > $PG/out/catchup_movies.log 2>&1 && rsync -a $PG/out/publish/movies/ $DST/movies/ --exclude "*_3d*"
    echo "$(date +%T) catch-up done rc=$? dumps $(ls data/red_*.npz | wc -l)"; } >> $CL 2>&1
  exit 0; fi
# pagefast-1102: the text rebuild (RSG pageH sync + mk_figs + build_page) runs as a small CPU step in the page's held allocation
# (pages_alloc.cfg) when called from the login node (USER: keep the login node free); inside a job/step it runs in place.
text_rebuild() {
  [ "$N" = rsgE ] && bash $J/rsg_E_1020/pageH/sync.sh   # 10-29 merged prodG+prodH view (pagefast-1102: ~6 s, was ~230 s)
  Q=$PG/out_quick; mkdir -p $Q/img $Q/publish; ln -sfn $DST/movies $Q/publish/movies; cd $PG/out || return 1
  ( ulimit -v 7000000; $PY $SRC/mk_figs.py data $Q/img > $Q/mk_figs.log 2>&1 ) || { tail -5 $Q/mk_figs.log; return 1; }
  cp -f img/movies.json img/fig_eq_latest.jpg $Q/img/ 2>/dev/null
  $PY $SRC/build_page.py $Q/img $Q/publish/index.html "page built $(date '+%Y-%m-%d %H:%M %Z') by quickpage.sh from the run's hst, logs and the $(ls data | grep -c '^red_.*\.npz$') dumps reduced so far" > $Q/build.log 2>&1 || { tail -5 $Q/build.log; return 1; }
  cp $Q/publish/index.html $DST/index.html
}
if [ "$2" = text ]; then text_rebuild; exit $?; fi
TX=
if [ -z "$SLURM_JOB_ID" ]; then
  while read pg jid node; do [ "$pg" = $N ] || continue
    [ "$(squeue -h -j $jid -o %T 2>/dev/null)" = RUNNING ] || continue
    timeout 900 srun --overlap --jobid=$jid -w $node -N1 -n1 -c 4 --mem=40G -J qp_$N nice -n 19 bash $J/quickpage.sh $N text < /dev/null && TX=ok || TX=fail
    break
  done < <(grep -v '^#' $J/pages_alloc.cfg)
fi
if [ "$TX" = fail ]; then exit 1; fi
if [ -z "$TX" ]; then text_rebuild || exit 1; fi
cd $PG/out || exit 1
# 10-24 USER "make it fast": the catch-up runs as an expansion job (catchup_job.sh, 24-wide), not 4-wide on login
NN=$(new | wc -l)
# 11-01 USER "make page updates fast and easy": catch-up runs DETACHED inside a held allocation (pages_alloc.cfg), never queued
if [ $NN -gt 0 ] && [ "$2" != nocatch ] && flock -n $PG/out/.catchup.lock true; then
  rm -f $DST/READY
  while read pg jid node; do [ "$pg" = $N ] || continue
    [ "$(squeue -h -j $jid -o %T 2>/dev/null)" = RUNNING ] || continue
    NSH=4 setsid nohup srun --overlap --jobid=$jid -w $node -N1 -n1 -c 4 --mem=60G nice -n 19 bash $J/catchup_job.sh $N > $PG/out/catchup_job.out 2>&1 < /dev/null &
    CA="in-alloc $jid $node"; break
  done < <(grep -v '^#' $J/pages_alloc.cfg)
  [ -z "$CA" ] && { sbatch -J pgcatch_$N $J/catchup_job.sh $N > /dev/null; CA="queued (no held alloc running)"; }
fi
echo "OK $N $(grep -oE 't = [0-9.e+]+ s = <b>[0-9.]+ d' $DST/index.html | sed 's/<b>//' | head -1) | $(du -k $DST/index.html | cut -f1) kB | missing figs $(grep -c 'missing figure' $DST/index.html) | $NN new dumps -> catch-up ${CA:-none}"
