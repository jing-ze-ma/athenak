#!/bin/bash
# accsp_prod_page job body (runs INSIDE an allocation via refresh.sh: srun --overlap, CPUs only, 8 cores, nice 10). READ ONLY on the run.
P=/resnick/groups/carnegie_poc/jingze/accsp_prod_page; J=/resnick/groups/carnegie_poc/jingze
export OUT=$P STATIC=${STATIC:-$J/ringsmooth_1102/long/L_gz} STREAM=${STREAM:-$J/accsp_1102/prod/fresh/stream}
source /etc/profile.d/modules.sh 2>/dev/null
module load python/3.11.6-gcc-11.3.1-rphh4kv py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui >/dev/null 2>&1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
FF=/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.3.1/ffmpeg-6.0-eim3de2umcr2h3ugaxw7as4pplxya727/bin/ffmpeg
NW=${NW:-4}   # pagefast-1102: 4 dumps x 2 threads = the 8 step cores (was 3); new frames render 4-wide in mk_figs (FRAME_PROCS)
T0=$SECONDS
cd $P; mkdir -p data img movie publish/img publish/movie
[ -s data/cs128_keband.json ] || python3 src/keband_cs128.py data/cs128_keband.json
echo "== reduce ($(python3 src/plan.py $STATIC/bin $STREAM/bin data | wc -l) dumps)"
python3 src/plan.py $STATIC/bin $STREAM/bin data | xargs -r -P $NW -L 1 python3 src/reduce_sp.py data
echo "== figs ($(( SECONDS - T0 )) s)"
python3 src/mk_figs.py || exit 1
python3 src/mk_trc.py || echo "mk_trc FAILED (page builds without the tracer section)"   # trc-1102
mkmp4() {  # <framedir> <name>
  local n fps; n=$(ls movie/$1/f_*.png 2>/dev/null | wc -l); [ "$n" -gt 0 ] || return 0
  if [ -f publish/movie/$2_frames.json ] && [ -z "$FORCE" ] && grep -q "\"n\": $n," publish/movie/$2_frames.json && \
     [ -z "$(find movie/$1 -name 'f_*.png' -newer publish/movie/$2.mp4)" ]; then return 0; fi
  fps=$(python3 -c "print(max(2, min(10, $n // 4)))")
  $FF -y -loglevel error -framerate $fps -pattern_type glob -i "movie/$1/f_*.png" -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,format=yuv420p" \
      -c:v libx264 -crf 23 -preset slow -movflags +faststart publish/movie/$2.tmp.mp4 && mv publish/movie/$2.tmp.mp4 publish/movie/$2.mp4
  cp "$(ls movie/$1/f_*.png | sort | tail -1)" publish/movie/$2_poster.png
  echo "{\"n\": $n, \"fps\": $fps}" > publish/movie/$2_frames.json
}
echo "== movies ($(( SECONDS - T0 )) s)"
mkmp4 fr_eq sp_eq; mkmp4 fr_ph sp_phot
echo "== page ($(( SECONDS - T0 )) s)"
python3 src/build_page.py || exit 1
echo "== done ($(( SECONDS - T0 )) s)"
