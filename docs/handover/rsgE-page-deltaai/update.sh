#!/bin/bash
# rsgE page (DeltaAI 10-20): one-command, incremental rebuild of the "RSG E DeltaAI" results page. READ ONLY on the run.
#   ./update.sh            -> OUT/publish/index.html + OUT/publish/movies/*.mp4 (publish the whole publish/ dir)
# Only new dumps are reduced (~45-60 s each incl. the T_gas EOS inversion, ~3.5 GB RSS, ~16 MB npz) and only new movie
# frames are rendered (4 movies, ~3 s per dump). Runs on the login node (nice 10).
# Override any of: RUN BINDIR EXTRA_RUN LOGDIR OUT ICF RPH0 JOBNAME SITE LAYOUT PY VISPY FFMPEG FPS CRF AGE MOVIE_RMAX
#   RUN defaults to /projects/.../run once its rsg.user.hst exists, else the /work/hdd run dir.
set -e
HERE=/work/nvme/bivj/jma20/rsgE_page
export SRC=$HERE/src
# Disk plan (10-20): link 1 writes /work/hdd/.../rsg_E_1021/run; afterwards a mover copies the dumps to /projects/.../run/bin and
# deletes the /work/hdd copies; link 2+ write /projects/.../run directly (hst copied there). Defaults follow that:
PROJ=/projects/bivj/jma20/rsg_E_1021/run
WORK=/work/hdd/bivj/jma20/rsg_E_1021/run
if [ -z "$RUN" ]; then
  if [ -f $PROJ/rsg.user.hst ] || [ ! -f $WORK/rsg.user.hst ]; then RUN=$PROJ; else RUN=$WORK; fi
fi
export RUN
export BINDIR=${BINDIR:-$PROJ/bin:$WORK/bin}          # ':'-list, first dir holding a dump index wins; missing dirs are fine
export EXTRA_RUN=${EXTRA_RUN:-$PROJ:$WORK}            # extra dirs searched for run.*.log (link 1 log stays in $WORK)
export LOGDIR=${LOGDIR:-/work/hdd/bivj/jma20/rsg_E_1021/logs:/projects/bivj/jma20/rsg_E_1021/logs}   # ':'-list of link.*.out dirs
export OUT=${OUT:-$HERE/out}
export ICF=${ICF:-/work/hdd/bivj/jma20/rsg_E_1021/files/rsg-E-1020/ic_rsg_g25_deep.txt}
[ -f "$ICF" ] || export ICF=$HERE/ic/ic_rsg_g25_deep.txt
[ -f "$ICF" ] || { mkdir -p $HERE/ic; git -C /u/jma20/ATHENAK/athenak show origin/bsg-files-1009:docs/handover/rsg-E-1020/ic_rsg_g25_deep.txt > $HERE/ic/ic_rsg_g25_deep.txt; }
export OLDRUN=""
export RPH0=${RPH0:-5.80019e13} JOBNAME=${JOBNAME:-rsgE21} AGE=${AGE:-120}
export VISPY=${VISPY:-/u/jma20/ATHENAK/builds/wt/hesnf_gpu/vis/python}
export PY=${PY:-/u/jma20/ATHENAK/venv/bin/python}
export SITE=${SITE:-DeltaAI}
export LAYOUT=${LAYOUT:-"2 nodes × 4 GH200, 24 MeshBlocks 548×64×64 (3 per GPU); 2 h ghx4-interactive links + 48 h ghx4 links"}
export FFMPEG=${FFMPEG:-$(/work/nvme/bivj/jma20/vtkvenv/bin/python -c "import imageio_ffmpeg;print(imageio_ffmpeg.get_ffmpeg_exe())")}
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MPLBACKEND=Agg
NICE="nice -n 10"
mkdir -p $OUT/data $OUT/img $OUT/publish/movies $OUT/logs
cd $OUT
exec 9>$OUT/.update.lock; flock -n 9 || { echo "update.sh already running"; exit 0; }
LOG=$OUT/logs/update.$(date +%Y%m%d_%H%M%S).log
{
  echo "update.sh RUN=$RUN BINDIR=$BINDIR LOGDIR=$LOGDIR OUT=$OUT"
  $NICE $PY $SRC/reduce.py $RUN data
  $NICE $PY $SRC/mk_figs.py data img > mk_figs.log 2>&1 || { tail -20 mk_figs.log; exit 1; }
  $NICE $PY $SRC/mk_movies.py data $OUT
  STATUS="page built $(date '+%Y-%m-%d %H:%M %Z') ($SITE) by update.sh from the run's hst, logs and the $(ls data/red_*.npz | wc -l) dumps reduced then"
  $NICE $PY $SRC/build_page.py img publish/index.html "$STATUS"
  ls -la publish/index.html publish/movies/
} 2>&1 | tee $LOG
exit ${PIPESTATUS[0]}
