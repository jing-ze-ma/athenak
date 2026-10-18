#!/bin/bash
# HGPAGE (10-18): build the "He Giant 2.75M Production" page for the hg145 production run. READ ONLY on the run.
# Portable: everything via env (defaults = Raven). One command:
#   RUN=<dir with hegiant.user.hst + run*.log> BINDIR=<dir with hegiant.{hydro_w,m1}.NNNNN.bin> OUT=<work dir> \
#   ICF=<ic_hegiant_14p5M_mk.txt> VISPY=<athenak>/vis/python PY=<python3 w/ numpy+matplotlib> \
#   PVPY=<pvbatch | python3 with the vtk wheel> FF=<ffmpeg with libx264> SITE=DeltaAI LAYOUT="..." ./upd_page.sh
# Deps: numpy, matplotlib (PY); VTK >= 9 with off-screen EGL or OSMesa (ParaView pvbatch, or `pip install vtk` / vtk-osmesa) (PVPY);
#       ffmpeg with libx264 (FF); AthenaK vis/python/bin_convert.py (VISPY). CPU only (no GPU needed), ~4 s per dump to reduce,
#       ~1 s per frame to render; run on a CPU compute/interactive node, not a login node, for long runs.
# Raven: ./upd_page.sh submits itself to partition interactive (8 cores, 60 GB, 1:59 h) and waits; INLINE=1 runs in place.
# FORCE=1 re-renders every frame (the colour scale also forces a full re-render by itself when it moves by > 5 %).
# Steps (incremental): reduce.py (new dumps -> OUT/data/red_NNNNN.npz) -> render3d.py (new frames -> OUT/movie/frames) ->
# ffmpeg (OUT/publish/movie3d/hg145.mp4, H.264 crf 22, 22 fps) + poster (last frame) -> mk_figs.py (OUT/img) -> build_page.py
# (OUT/publish/index.html, figures inline). Publish index.html with files movie3d/hg145.mp4 + movie3d/hg145_poster.png.
SRC=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-/raven/ptmp/jinma/hg145_prod_1018/page}
if [ -z "$SLURM_JOB_ID" ] && [ -z "$INLINE" ] && [ -z "$NO_SBATCH" ] && [[ $(hostname) == raven* ]]; then
  mkdir -p $OUT/logs
  exec sbatch --wait -p ${PART:-interactive} -n 1 -c 8 --mem=60G -t 01:59:00 -J hgpage -o $OUT/logs/upd.%j.out --export=ALL "$SRC/upd_page.sh"
fi
set -e
export RUN=${RUN:-/raven/ptmp/jinma/hg145_prod_1018/prod/R}
export BINDIR=${BINDIR:-$RUN/bin}
export ICF=${ICF:-/raven/ptmp/jinma/ic_star_runs/hegiant_14p5M_mk/ic_hegiant_14p5M_mk.txt}
export VISPY=${VISPY:-/raven/u/jinma/ATHENAK/wt-accbig/vis/python}
export SITE=${SITE:-Raven} LAYOUT=${LAYOUT:-"1 node × 4 A100, 4 MeshBlocks 704×64×64; chain of 24 h links"}
PY=${PY:-/u/jinma/miniconda3/envs/magritte/bin/python}
FF=${FF:-/u/jinma/miniconda3/envs/magritte/bin/ffmpeg}
PVPY=${PVPY:-/mpcdf/soft/SLE_15/packages/x86_64/paraview/6.1.1/bin/pvbatch}
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MPLBACKEND=Agg
mkdir -p $OUT/data $OUT/img $OUT/movie/frames $OUT/publish/movie3d
cd $OUT
[ -f publish/movie3d/hg145.json ] || cp $SRC/hg145.json publish/movie3d/hg145.json
$PY $SRC/reduce.py $RUN data
$PVPY $SRC/render3d.py data movie/frames publish/movie3d/hg145.json 2>&1 | grep -vE "bad X server|Deprecat|AddActor2D" || true
cp movie/frames/scale.json publish/movie3d/hg145_scale.json
FPS=$($PY -c "import json;print(json.load(open('publish/movie3d/hg145.json'))['fps'])")
CRF=$($PY -c "import json;print(json.load(open('publish/movie3d/hg145.json'))['crf'])")
ls $OUT/movie/frames/f_*.png | sort > movie/list.txt
NF=$(wc -l < movie/list.txt)
$FF -y -loglevel error -framerate $FPS -pattern_type glob -i "$OUT/movie/frames/f_*.png" -vf "format=yuv420p" -c:v libx264 \
    -crf $CRF -preset slow -movflags +faststart publish/movie3d/hg145.tmp.mp4
mv publish/movie3d/hg145.tmp.mp4 publish/movie3d/hg145.mp4
cp "$(tail -1 movie/list.txt)" publish/movie3d/hg145_poster.png
echo "{\"n\": $NF}" > publish/movie3d/hg145_frames.json
$PY $SRC/mk_figs.py data img > mk_figs.log 2>&1 || { tail -20 mk_figs.log; exit 1; }
STATUS="page built $(date '+%Y-%m-%d %H:%M %Z') ($SITE) by upd_page.sh from the run's hst, logs and the $NF dumps present then"
$PY $SRC/build_page.py img publish/index.html "$STATUS"
ls -la publish/index.html publish/movie3d/
