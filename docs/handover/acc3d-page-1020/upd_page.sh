#!/bin/bash
# ACC3D page bundle 10-20 (portable): build the "RY Per Accretor 3-D" page from a 3-D accretor run. READ ONLY on the run.
# Everything via env (defaults = Raven run12). One command, e.g. Caltech:
#   RUN=/resnick/groups/carnegie_poc/jingze/acc3d_trc_1019 OUT=<work dir> PY=<python3 with numpy+matplotlib> \
#   FF=<ffmpeg with libx264> PVPY=<pvbatch | python3 with vtk>=9 (offscreen)> EOSTBL=<eos table dump> \
#   SITE=Caltech LAYOUT="1 node x 2 H200, 8 MeshBlocks 420x58x64 (4 per GPU)" BIN_LABEL="..." ./upd_page.sh
# Env:
#   RUN      run parent with static/ and stream/ (ryper.{user,hydro,trc}.hst, run.log, bin/ryper.{hydro_w,m1}.NNNNN.bin, rst/)
#   OUT      work dir (data caches, img/, movie frames; publish/ = what gets published)
#   PY FF PVPY  python3 (numpy, matplotlib), ffmpeg (libx264), pvbatch or a python3 with vtk >= 9 off-screen (EGL/OSMesa)
#   VISPY    optional: an AthenaK vis/python dir whose bin_convert.read_binary replaces the bundled binread.py
#   EOSTBL   the run's EOS table dump (eos_x0.70_y0.28_logd-18.txt, md5 e548d9c2...; mk_eostbl.sh regenerates it)
#   OPAC     Rosseland table of the deck (default: ../acc3d-prod-1018/rosseland_mesakap_gs98_x0.70_z0.02.txt in this checkout)
#   SITE LAYOUT BIN_LABEL DECK_LABEL RUN_LABEL CHAIN_LABEL   page text;  JOBNAME  slurm job name of the links (queue line)
#   NP (4) worker processes; NTRC (4) tracer-analysis times; FORCE=1 redraw everything; SKIP3D=1 / SKIP2D=1 skip a movie
#   Raven only: without INLINE/NO_SBATCH/SLURM_JOB_ID the script submits itself to partition interactive (8 cores, 60 GB, 1:59 h).
# Steps (incremental, caches in $OUT/data): mk_prod.py (status, checks, budgets, history, maps) -> mk_tracer.py (trc.hst + s0 analyses)
#   -> mk_omega.py -> mk_movie.py (2-D frames) -> reduce3d.py + render3d.py (3-D frames) -> ffmpeg -> build_page.py.
# Output: $OUT/publish/index.html + publish/movie/acc3d_{eq,3d}.mp4 + acc3d_{eq,3d}_poster.png (publish these, same paths).
export SRC=${SRC:-$(cd "$(dirname "$0")" && pwd)}
export OUT=${OUT:-/raven/ptmp/jinma/acc3d_1017/page12}
if [ -z "$SLURM_JOB_ID" ] && [ -z "$INLINE" ] && [ -z "$NO_SBATCH" ] && [[ $(hostname) == raven* ]]; then
  mkdir -p $OUT/logs
  exec sbatch --wait -p ${PART:-interactive} -n 1 -c 8 --mem=60G -t 01:59:00 -J accpage -o $OUT/logs/upd.%j.out --export=ALL "$SRC/upd_page.sh"
fi
set -e
export RUN=${RUN:-/raven/ptmp/jinma/acc3d_1017/prod/run12}
export EOSTBL=${EOSTBL:-/raven/ptmp/jinma/acc3d_1017/mk/eosdump/eos_x0.70_y0.28_logd-18.txt}
export OPAC=${OPAC:-$SRC/../acc3d-prod-1018/rosseland_mesakap_gs98_x0.70_z0.02.txt}
export SITE=${SITE:-Raven} LAYOUT=${LAYOUT:-"2 nodes x 4 A100, 8 MeshBlocks 420x58x64; 24 h links"}
export BIN_LABEL=${BIN_LABEL:-"prodmerge-1019b @ 13206541 (A100 nofma md5 29caa1ab)"} DECK_LABEL=${DECK_LABEL:-acc3d_prod_v11amwtr}
export RUN_LABEL=${RUN_LABEL:-"fresh tracer run (run12)"} CHAIN_LABEL=${CHAIN_LABEL:-"prod_chain links, fresh from t = 0"}
export NP=${NP:-4} JOBNAME=${JOBNAME:-acc3dP}
PY=${PY:-/u/jinma/miniconda3/envs/magritte/bin/python}
FF=${FF:-/u/jinma/miniconda3/envs/magritte/bin/ffmpeg}
PVPY=${PVPY:-/mpcdf/soft/SLE_15/packages/x86_64/paraview/6.1.1/bin/pvbatch}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
[ -d "$RUN/static" ] || [ -d "$RUN/stream" ] || { echo "RUN=$RUN has no static/ or stream/"; exit 1; }
mkdir -p $OUT/data $OUT/img $OUT/logs $OUT/movie/frames2d $OUT/movie/frames3d $OUT/publish/movie
cd $OUT
T0=$SECONDS
step() { echo "== $1 ($(( SECONDS - T0 )) s)"; }
step mk_prod;   $PY $SRC/mk_prod.py img   > logs/mk_prod.log 2>&1   || { tail -30 logs/mk_prod.log; exit 1; }
step mk_tracer; $PY $SRC/mk_tracer.py img > logs/mk_tracer.log 2>&1 || { tail -30 logs/mk_tracer.log; exit 1; }
step mk_omega;  $PY $SRC/mk_omega.py img  > logs/mk_omega.log 2>&1  || { tail -30 logs/mk_omega.log; exit 1; }
mkmp4() {  # <framedir> <name> <fps> <crf>
  local n; n=$(ls $1/f_*.png 2>/dev/null | wc -l); [ "$n" -gt 0 ] || return 0
  $FF -y -loglevel error -framerate $3 -pattern_type glob -i "$1/f_*.png" -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,format=yuv420p" \
      -c:v libx264 -crf $4 -preset slow -movflags +faststart publish/movie/$2.tmp.mp4
  mv publish/movie/$2.tmp.mp4 publish/movie/$2.mp4
  cp "$(ls $1/f_*.png | sort | tail -1)" publish/movie/$2_poster.png
  echo "{\"n\": $n}" > publish/movie/$2_frames.json
}
if [ -z "$SKIP2D" ]; then
  step mk_movie; $PY $SRC/mk_movie.py movie/frames2d > logs/mk_movie.log 2>&1 || { tail -30 logs/mk_movie.log; exit 1; }
  tail -1 logs/mk_movie.log; mkmp4 movie/frames2d acc3d_eq 10 23
fi
if [ -z "$SKIP3D" ]; then
  step reduce3d; $PY $SRC/reduce3d.py data/red3d > logs/reduce3d.log 2>&1 || { tail -30 logs/reduce3d.log; exit 1; }
  tail -1 logs/reduce3d.log
  step render3d
  for w in $(seq 0 $((NP-1))); do
    $PVPY $SRC/render3d.py data/red3d movie/frames3d $SRC/mov3d.json $w/$NP > logs/render3d.$w.log 2>&1 &
  done
  wait
  grep -hE "Error|error" logs/render3d.*.log | grep -v "bad X server" | head -5 || true
  FPS=$($PY -c "import json;print(json.load(open('$SRC/mov3d.json'))['fps'])"); CRF=$($PY -c "import json;print(json.load(open('$SRC/mov3d.json'))['crf'])")
  mkmp4 movie/frames3d acc3d_3d $FPS $CRF
fi
step build_page
STATUS="page built $(date '+%Y-%m-%d %H:%M %Z') on $SITE by acc3d-page-1020/upd_page.sh from the hst, logs and dumps present then"
$PY $SRC/build_page.py img publish/index.html "$STATUS"
step done
ls -la publish/index.html publish/movie/
