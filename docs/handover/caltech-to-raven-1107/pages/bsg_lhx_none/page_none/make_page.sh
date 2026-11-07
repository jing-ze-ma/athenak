#!/bin/bash
# bsgli_1104 page (none): Linhao-exact BSG run /resnick/groups/carnegie_poc/jingze/bsgli_1104/prod (linhao-prod-1104, binary 26eae921), copy of
# bsgli_1014/pagework_li.  Page rule USER 10-24: no inlining; publish index.html + media files.
# usage: make_page.sh [pagedir]   CPU only; run inside a held allocation (refresh.sh does that).
export W=/resnick/groups/carnegie_poc/jingze/bsgli_1104/page_none
PAGE=${1:-/resnick/home/jingze/ATHENAK/pages/bsg_lhx_none}
source /etc/profile >/dev/null 2>&1; module purge
module load python/3.11.6-gcc-11.3.1-rphh4kv py-scipy/1.11.1-gcc-11.3.1-swvnmdi \
  py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui >/dev/null 2>&1
source /resnick/groups/carnegie_poc/jingze/bsg_1009/pagework/venv/bin/activate            # system site + astropy (Fig 6)
export PATH=$PATH:/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.3.1/ffmpeg-6.0-eim3de2umcr2h3ugaxw7as4pplxya727/bin
export FFMPEG_MODULE=none
export BSG_REPRO_RUN=/resnick/groups/carnegie_poc/jingze/bsgli_1104/prod
export BSG_REPRO_FULL_RUN=$BSG_REPRO_RUN BSG_BUILD_MERGED=
export BSG_INP=/resnick/groups/carnegie_poc/jingze/bsgli_1104/bsgLI_exact.athinput
# IC overlays (Fig 3 right panel, r_photo,IC, Fe peak): L. Ma's profile (the run's he_ic_linhao_file)
# converted to the 7-column layout + his isothermal extension above rmax (tools/make_ic_linhao.py)
export BSG_IC=$W/ic/ic_linhao_7col.txt
export BSG_OPAC=/resnick/groups/carnegie_poc/jingze/bsgli_1014/files/rosseland_ma2026_x0.7_z0.008_padded.txt
export BSG_COST_PARTS=$W/parts.json BSG_RUN_LABEL=bsgLX-none
export BSG_KE_LABEL='Linhao-exact (linhao-prod-1104), inner-face dissipation none'
export BSG_HW_LABEL='1 Caltech node: 4 H200 (4 MeshBlocks per GPU), hpc-sm-01-05'
export MPLCONFIGDIR=$W/pipe/cache/mpl
mkdir -p $PAGE/img $PAGE/movie $MPLCONFIGDIR $W/ic
[ -s $BSG_IC ] || python3 $W/tools/make_ic_linhao.py /resnick/groups/carnegie_poc/jingze/bsgli_1014/files/linhao_profile_full.txt $BSG_IC || exit 1
python3 -c "import astropy" || echo "WARNING: astropy missing: Fig 6 / bsg_figs will fail"
# 57.4-d comparison: reduce every dump with 30 <= t <= 57.41 d not yet in cache/cmp57
NS=$(python3 - <<'PY'
import glob, os, re
run = os.environ['BSG_REPRO_RUN']
out = []
for h in sorted(glob.glob(run + '/bin/bsg3d.hydro_w.*.bin')):
    n = int(h.split('.')[-2])
    t = float(re.search(r'time=\s*([-+0-9.eE]+)', open(h, 'rb').read(400).decode('latin-1')).group(1)) / 86400
    if 30 <= t <= 57.41 and not os.path.exists('%s/pipe/cache/cmp57/c_%05d.npz' % (os.environ['W'], n)):
        out.append(str(n))
print(' '.join(out))
PY
)
[ -n "$NS" ] && { echo "cmp57_dumps $NS"; python3 $W/pipe/cmp57_dumps.py $NS || echo "cmp57_dumps FAILED"; }
cp $W/pipe/index.html $PAGE/index.html
bash $W/pipe/update_page_repro.sh $PAGE
python3 $W/tools/dark_paper_fig4.py $PAGE   # 10-14: img/paper_fig4_dark.png (Fig 4 dark panel)
python3 $W/tools/inject_status.py $PAGE
