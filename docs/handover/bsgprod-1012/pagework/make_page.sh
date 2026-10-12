#!/bin/bash
# Regenerate the BSG true-reproduction page (Caltech run bsg_1009/prod_noseed, read only) with viper's
# page pipeline (fork branch bsg-files-1009, copy in ./pipe with Caltech edits; viper's original
# index.html kept as pipe/index.viper.html).  Run on a CPU compute node, e.g.
#   sbatch -A carnegie_poc -p expansion -c 8 --mem=64G -t 02:00:00 -o make_page.%j.out make_page.sh
# usage: make_page.sh [pagedir]   (default /resnick/home/jingze/ATHENAK/pages/bsg_truerepro)
# Writes only under ./pipe (out/ cache/ movie/) and pagedir.  Paper crops img/paper_fig*.png are
# not made here (copied from the old artifact at publish time).
export W=/resnick/groups/carnegie_poc/jingze/bsg_1009/pagework
PAGE=${1:-/resnick/home/jingze/ATHENAK/pages/bsg_truerepro}
source /etc/profile >/dev/null 2>&1; module purge
module load python/3.11.6-gcc-11.3.1-rphh4kv py-scipy/1.11.1-gcc-11.3.1-swvnmdi \
  py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui >/dev/null 2>&1
source $W/venv/bin/activate            # system site + astropy (Fig 6)
export PATH=$PATH:/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.3.1/ffmpeg-6.0-eim3de2umcr2h3ugaxw7as4pplxya727/bin
export FFMPEG_MODULE=none
export BSG_REPRO_RUN=/resnick/groups/carnegie_poc/jingze/bsg_1009/prod_noseed
export BSG_REPRO_FULL_RUN=$BSG_REPRO_RUN BSG_BUILD_MERGED=
export BSG_INP=/resnick/groups/carnegie_poc/jingze/bsg_1009/files/bsg3d_truerepro2_hr_lm.athinput
export BSG_IC=/resnick/groups/carnegie_poc/jingze/bsg_1009/files/ic_ma2026_f0.txt
export BSG_OPAC=/resnick/groups/carnegie_poc/jingze/bsg_1009/files/rosseland_ma2026_x0.7_z0.008_padded.txt
export BSG_COST_PARTS=$W/parts.json BSG_RUN_LABEL=truerepro2
export BSG_HW_LABEL='1 Caltech node, 2 H200, 16 MeshBlocks (8 per GPU)'
export MPLCONFIGDIR=$W/pipe/cache/mpl
mkdir -p $PAGE/img $PAGE/movie $MPLCONFIGDIR
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
