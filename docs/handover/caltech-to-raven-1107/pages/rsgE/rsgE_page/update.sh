#!/bin/bash
# rsgE page (Caltech copy of the DeltaAI 10-20 pipeline, 10-21): one-command, incremental rebuild of the "RSG E DeltaAI" results page. READ ONLY on the run.
#   ./update.sh            -> OUT/publish/index.html + OUT/publish/movies/*.mp4 (publish the whole publish/ dir)
# Only new dumps are reduced (~45-60 s each incl. the T_gas EOS inversion, ~3.5 GB RSS, ~16 MB npz) and only new movie
# frames are rendered (4 movies, ~3 s per dump). Run it as a CPU job (expansion), not on the login node.
# Override any of: RUN BINDIR EXTRA_RUN LOGDIR OUT ICF RPH0 JOBNAME SITE LAYOUT PY VISPY FFMPEG FPS CRF AGE MOVIE_RMAX
#   RUN defaults to /projects/.../run once its rsg.user.hst exists, else the /work/hdd run dir.
set -e
HERE=/resnick/groups/carnegie_poc/jingze/rsgE_page
export SRC=$HERE/src
R=/resnick/groups/carnegie_poc/jingze/rsg_E_1020
export RUN=${RUN:-$R/pageH} BINDIR=${BINDIR:-$R/pageH/bin} EXTRA_RUN=${EXTRA_RUN:-$R/prodH} LOGDIR=${LOGDIR:-$R/prodH}   # 10-29: pageH = prodG < 1676 d + prodH (re-run from 1675 d with vtxfix); prodG >= 1676 d caches in out/*_prodG_post1675
export OUT=${OUT:-$HERE/out} ICF=${ICF:-$HERE/ic/ic_rsg_g25_deep.txt}
export OLDRUN=""
export RPH0=${RPH0:-5.80019e13} JOBNAME=${JOBNAME:-rsgED,rsgEB,rsgEH,rsgEH1} AGE=${AGE:-120}
export VISPY=${VISPY:-/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d}
export PY=${PY:-/resnick/groups/carnegie_poc/jingze/venv_vtk/bin/python}   # + scipy, h5py (pip 10-21)
export SITE=${SITE:-Caltech}
export LAYOUT=${LAYOUT:-"1 node × 4 H200: 30-min debug-QOS links until 75 d, then an 8-h backfill job (2-h chunks, a restart each) until a longer allocation starts; 24 MeshBlocks 548×64×64 (6 per GPU)"}
export FFMPEG=${FFMPEG:-/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.3.1/ffmpeg-6.0-eim3de2umcr2h3ugaxw7as4pplxya727/bin/ffmpeg}
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MPLBACKEND=Agg
NICE="nice -n 10"
mkdir -p $OUT/data $OUT/img $OUT/publish/movies $OUT/logs
cd $OUT
exec 9>$OUT/.update.lock; flock -n 9 || { echo "update.sh already running"; exit 3; }   # 10-22: nonzero so callers do not stage a stale publish/
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
