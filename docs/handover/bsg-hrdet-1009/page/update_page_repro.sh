#!/bin/bash
# One-shot refresh of the BSG REPRODUCTION page (run $BSG_REPRO_RUN, default repro_sf_4n since
# 10-06; the unseeded full run $BSG_REPRO_FULL_RUN = repro_4n only for the KE panel; read only):
#   figures (bsg_figs.py) -> cost (cost_repro.py) -> status.json (gen_status.py) -> movie frames
#   + mp4 (movie/make_movie.sh) -> copy into <pagedir> -> <pagedir>/new_files.txt (what to publish).
# Colour ranges are re-chosen from the newest dump on every refresh and all frames redrawn with
# them (BSG_RECHOOSE=1, fixed within each movie); set BSG_RECHOOSE=0 to keep ranges.json and
# render only new frames.  Login node: nice, 8 threads.  Writes only under scripts_repro/ and
# <pagedir>.   usage: update_page_repro.sh [pagedir]
S=$(cd "$(dirname "$0")" && pwd)   # handover copy: the script's own dir
SP=${1:-/viper/ptmp2/jinma/pages_1004/bsgrepro}
export BSG_RECHOOSE=${BSG_RECHOOSE:-1}
# BSG_REPRO_RUN before 10-07: /raven/ptmp/jinma/bsg_raven_1002/repro_sf_4n (Raven only); since
# 10-07 the merged read-only view Raven (t <= t_rst 26.49 d) + viper (bsg_viper_1006/repro),
# rebuilt below by build_merged.py (cache/merged_repro -> cache/repro_sf_4n symlink).
export BSG_REPRO_RUN=${BSG_REPRO_RUN:-/viper/ptmp2/jinma/bsg_viper_1006/merged_repro}
export BSG_REPRO_FULL_RUN=${BSG_REPRO_FULL_RUN:-/raven/ptmp/jinma/bsg_raven_1002/repro_4n}
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd $S || exit 1
mkdir -p $SP/img $SP/movie cache out
[ -e cache/merged_repro ] || ln -s repro_sf_4n cache/merged_repro
# viper only (merged Raven+viper view); skipped when BSG_BUILD_MERGED is empty
BSG_BUILD_MERGED=${BSG_BUILD_MERGED-/viper/ptmp2/jinma/bsg_viper_1006/build_merged.py}
if [ -n "$BSG_BUILD_MERGED" ]; then
python3 $BSG_BUILD_MERGED repro > cache/merged_run.log 2>&1 || {
  echo "build_merged FAILED"; tail -3 cache/merged_run.log; exit 1; }
fi
nice -n 10 timeout 7000 python3 $S/bsg_figs.py --arm repro > cache/repro_run.log 2>&1 || {
  echo "bsg_figs FAILED"; tail -3 cache/repro_run.log; exit 1; }
nice -n 10 python3 $S/cost_repro.py > cache/cost_run.log 2>&1 || {
  echo "cost_repro FAILED"; tail -3 cache/cost_run.log; exit 1; }
nice -n 10 python3 $S/ke_compare.py > cache/ke_run.log 2>&1 || {
  echo "ke_compare FAILED"; tail -3 cache/ke_run.log; exit 1; }
for f in ours_fig2 ours_fig3 ours_fig4 ours_fig6 ours_int ours_ke; do cp out/$f.png $SP/img/; done
python3 $S/gen_status.py $SP > /dev/null || { echo "gen_status FAILED"; exit 1; }
# 10-08: the 57.4 d comparison section (needs cache/cmp57 from cmp57_dumps.py)
nice -n 10 python3 $S/compare57.py $SP > $S/cache/cmp57_agg.log 2>&1 || echo "compare57 FAILED"
cp $S/out/ours_cmp57.png $SP/img/
nice bash $S/movie/make_movie.sh > $S/movie/last_make.log 2>&1 || echo "make_movie FAILED (see $S/movie/last_make.log)"
# new_files.txt: every file of the page that changed in this refresh (publish list)
: > $SP/new_files.txt
for f in index.html status.json img/ours_fig2.png img/ours_fig3.png img/ours_fig4.png \
         img/ours_fig6.png img/ours_int.png img/ours_ke.png img/ours_cmp57.png; do echo $f >> $SP/new_files.txt; done
for k in f3 f4 fint; do
  mkdir -p $SP/movie/$k
  for f in $S/movie/frames/$k/*.png; do
    b=$(basename $f)
    if ! cmp -s $f $SP/movie/$k/$b; then cp $f $SP/movie/$k/; echo movie/$k/$b >> $SP/new_files.txt; fi
  done
done
for m in fig3 fig4 figint; do
  [ -f $S/movie/$m.mp4 ] && ! cmp -s $S/movie/$m.mp4 $SP/movie/$m.mp4 && {
    cp $S/movie/$m.mp4 $SP/movie/; echo movie/$m.mp4 >> $SP/new_files.txt; }
done
N=$(ls $SP/movie/f4 | wc -l)
# n.json: frame count + each frame's dump time [d] (dumps are not exactly 0.5 d apart after
# the restarts; the page's movie label reads t[i], 10-08)
python3 $S/movie_times.py "$SP" "$BSG_REPRO_RUN"
echo movie/n.json >> $SP/new_files.txt
for p in paper_fig2 paper_fig3 paper_fig4 paper_fig6; do
  [ -f $SP/img/$p.png ] && echo img/$p.png >> $SP/new_files.txt; done
python3 - "$SP" <<'PY'
import json, re, sys
d = json.load(open(sys.argv[1] + '/status.json'))
print('OK |', d['s-time'], '| wall', d['c-wall'], 'proj', d['c-proj'], '| fig4',
      re.sub('<[^>]+>', '', d['w-fig4']).replace('&#9737;', 'sun'))
PY
echo "movie frames: $N (files listed in $SP/new_files.txt: $(wc -l < $SP/new_files.txt))"
