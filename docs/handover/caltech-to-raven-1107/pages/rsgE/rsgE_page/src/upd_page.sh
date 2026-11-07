#!/bin/bash
# RSGPAGE (10-19): build the "RSG g25 Production" page for the rsg_g25 cubed-sphere production run. READ ONLY on the run.
# Env (defaults = Raven): RUN (dir with rsg.user.hst, run.log, link.*.out), BINDIR (rsg.{hydro_w,m1}.NNNNN.bin), OUT (work dir),
#   ICF (IC column), OLDRUN (archived old-seed run, for the KE comparison), VISPY (AthenaK vis/python), PY (numpy+matplotlib).
# Raven: ./upd_page.sh submits itself to partition interactive (8 cores, 60 GB) and waits; INLINE=1 runs in place.
# Steps (incremental): reduce.py (new 1 d dumps -> OUT/data/red_NNNNN.npz, ~7 s per dump) -> mk_figs.py (OUT/img)
#   -> build_page.py (OUT/publish/index.html, figures inlined base64, single file). Publish index.html alone.
export SRC=${SRC:-$(cd "$(dirname "$0")" && pwd)}
OUT=${OUT:-/raven/ptmp/jinma/rsg_g25_1019/page}
if [ -z "$SLURM_JOB_ID" ] && [ -z "$INLINE" ] && [ -z "$NO_SBATCH" ] && [[ $(hostname) == raven* ]]; then
  mkdir -p $OUT/logs
  exec sbatch --wait -p ${PART:-interactive} -n 1 -c 8 --mem=60G -t 01:00:00 -J rsgpage -o $OUT/logs/upd.%j.out --export=ALL "$SRC/upd_page.sh"
fi
set -e
export RUN=${RUN:-/raven/ptmp/jinma/rsg_g25_1019/prod/R}
export BINDIR=${BINDIR:-$RUN/bin}
export ICF=${ICF:-/raven/ptmp/jinma/rsg_g25_1019/ic/pass3b/ic_rsg_g25.txt}
export OLDRUN=${OLDRUN:-/raven/ptmp/jinma/rsg_g25_1019/prod/R_oldseed}
export VISPY=${VISPY:-/raven/u/jinma/ATHENAK/wt-accbig/vis/python}
export SITE=${SITE:-Raven} LAYOUT=${LAYOUT:-"6 nodes × 4 A100, 24 MeshBlocks 591×64×64 (1 per GPU); chain of 24 h links"}
PY=${PY:-/u/jinma/miniconda3/envs/magritte/bin/python}
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MPLBACKEND=Agg
mkdir -p $OUT/data $OUT/img $OUT/publish
cd $OUT
$PY $SRC/reduce.py $RUN data
$PY $SRC/mk_figs.py data img > mk_figs.log 2>&1 || { tail -20 mk_figs.log; exit 1; }
STATUS="page built $(date '+%Y-%m-%d %H:%M %Z') ($SITE) by upd_page.sh from the run's hst, logs and the $(ls data/red_*.npz | wc -l) dumps reduced then"
$PY $SRC/build_page.py img publish/index.html "$STATUS"
ls -la publish/index.html
