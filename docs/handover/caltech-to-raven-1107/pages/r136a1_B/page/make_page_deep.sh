#!/bin/bash
# r136a1-1104 DEEP run page: figures + movie frames + status.json from runB/prod (new dumps only; --full redraws all frames
# with re-chosen colour ranges). CPU only, inside a held allocation (refresh.sh). Page rule USER 10-24: no inlining,
# index.html + media files.   usage: make_page.sh PAGEDIR [--full]
W=/resnick/groups/carnegie_poc/jingze/r136a1_1104/page; PAGE=${1:-/resnick/home/jingze/ATHENAK/pages/r136a1_deep}
export R136_RUN=/resnick/groups/carnegie_poc/jingze/r136a1_1104/runB/prod R136_GM0=3.1017784815e28 R136_IC=/resnick/groups/carnegie_poc/jingze/r136a1_1104/ic/ic_r136a1_B.txt R136_LABEL="R136a1 phase 1 deep B (128\$^2\$)" R136_TLIM=${R136_TLIM:-1.728e6} R136_BANNERHEAD="Phase 1, deep inner boundary (x1min 19.45 Rsun, log T 6.5), 128&times;128 lateral."
source /etc/profile >/dev/null 2>&1; module purge
module load python/3.11.6-gcc-11.3.1-rphh4kv py-scipy/1.11.1-gcc-11.3.1-swvnmdi py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui >/dev/null 2>&1
export MPLCONFIGDIR=$W/cache/mpl OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
mkdir -p $PAGE/img $PAGE/movie $MPLCONFIGDIR
cp $W/index_deep.html $PAGE/index.html
python3 $W/r136_page.py $PAGE $2
# 11-05 USER: ADDITIONAL HG/RSG-style line plots (img/hgstyle_*.png); the figures above are kept unchanged
python3 $W/hgstyle_figs.py $PAGE || echo "hgstyle_figs FAILED (page figures above are unaffected)"
