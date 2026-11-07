#!/bin/bash
# hgglob page (10-21): one-command, incremental rebuild of the "He Giant Global Caltech" results page for the He giant GLOBAL
# cubed-sphere run /resnick/groups/carnegie_poc/jingze/hgglob_1021 (READ ONLY on the run). Pipeline = copy of rsgE_page.
#   sbatch -A carnegie_poc -p expansion -c 8 --mem=32G -t 01:30:00 -J hgGpage -o <OUT>/logs/page.%j.out --wrap "bash update.sh"
#   -> OUT/publish/index.html + OUT/publish/movies/*.mp4 (publish the whole publish/ dir)
# Only new dumps are reduced (~1-1.5 min each incl. the T_gas EOS inversion, ~5 GB RSS) and only new movie frames are rendered.
# Run it as a CPU job (expansion), not on the login node.
# Override any of: RUN BINDIR EXTRA_RUN LOGDIR OUT ICF RPH0 BASENAME JOBNAME SITE LAYOUT PY VISPY FFMPEG FPS CRF AGE MOVIE_RMAX TFLOOR
set -e
HERE=/resnick/groups/carnegie_poc/jingze/hgglob_page
export SRC=$HERE/src
R=/resnick/groups/carnegie_poc/jingze/hgglob_1021
export RUN=${RUN:-$R/prodG} BINDIR=${BINDIR:-$R/prodG/bin} EXTRA_RUN=${EXTRA_RUN:-$R/prodG} LOGDIR=${LOGDIR:-$R/prodG}   # 10-22: FRESH N^2-seed run (old prodF -> out/old_prodF_1022)
export OUT=${OUT:-$HERE/out} ICF=${ICF:-$R/files/hg145-prod-1018/ic_hegiant_14p5M_mk.txt}
export OLDRUN=""
export BASENAME=${BASENAME:-hegiant}
export RPH0=${RPH0:-7.94653e12}            # IC column header R_ph (ic_hegiant_14p5M_mk.txt)
export RIN=${RIN:-5.5656e11}               # deck <mesh>/x1min (inner wall), drawn dashed in the zooms
export TFLOOR=${TFLOOR:-2400}              # deck <hydro>/tfloor_kelvin
export JOBNAME=${JOBNAME:-hgGD,hgGS,hgGH3,hgGH1,hgGB} AGE=${AGE:-120}
export VISPY=${VISPY:-/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d}
export PY=${PY:-/resnick/groups/carnegie_poc/jingze/venv_vtk/bin/python}   # numpy, scipy, matplotlib, h5py, PIL
export SITE=${SITE:-Caltech}
export LAYOUT=${LAYOUT:-"1 node × 4 H200 30-min debug-QOS links (6 MeshBlocks per GPU) until the 3-node 7-d allocation starts (est. 22:49); 24 MeshBlocks 704×64×64"}
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
