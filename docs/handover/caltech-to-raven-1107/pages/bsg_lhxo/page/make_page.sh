#!/bin/bash
# 11-05 ~19:10: REPOINTED to prod_vsnmlat_ns (vstar_n + merge_lat x2 below 35 Rsun, NO sink; deck bsgLI_exact_vsn_mlat md5 1126c821; binary doface_v2 bd820f7f;
# from t=0 ~20:15, tlim 57.4 d). The vsn+sink arm (prod_vsnsink, stopped ~1.2 d) is dropped; prod_doface (vstar_n only, ~18 d) stays a comparison curve.
# 11-05 evening: REPOINTED to prod_vsnsink (vstar_n + he_athpp_wall_sink_L 2.79e40, dtref 212; binary athena_gpu_he_star_m1_doface_v2_fmadoff md5 bd820f7f,
# from t=0 at 18:21, tlim 57.4 d). The stopped vstar_n-only arm prod_doface (~18 d) is now a comparison curve in pipe/cmp_arms.py.
# doface_1105 page: BSG Athena++-equivalent arm (Linhao-exact deck + rad_m1/doppler_vel = vstar_n) /resnick/groups/carnegie_poc/jingze/doface_1105/prod_doface
# (athpp-doface-1105, binary athena_gpu_he_star_m1_doface_v1_fmadoff md5 fcdbd608), clone of doface_1105/page (11-05).
# Extra vs page_none: pipe/cmp_arms.py (this arm vs the Linhao-exact 'none' production bsgli_1104/prod: r_photo, near-base Mdot/M_env, L, wall-cell T).  Page rule USER 10-24: no inlining; publish index.html + media files.
# usage: make_page.sh [pagedir]   CPU only; run inside a held allocation (refresh.sh does that).
export W=/resnick/groups/carnegie_poc/jingze/doface_1105/page
PAGE=${1:-/resnick/home/jingze/ATHENAK/pages/bsg_doface}
source /etc/profile >/dev/null 2>&1; module purge
module load python/3.11.6-gcc-11.3.1-rphh4kv py-scipy/1.11.1-gcc-11.3.1-swvnmdi \
  py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui >/dev/null 2>&1
source /resnick/groups/carnegie_poc/jingze/bsg_1009/pagework/venv/bin/activate            # system site + astropy (Fig 6)
export PATH=$PATH:/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.3.1/ffmpeg-6.0-eim3de2umcr2h3ugaxw7as4pplxya727/bin
export FFMPEG_MODULE=none
export BSG_REPRO_RUN=/resnick/groups/carnegie_poc/jingze/doface_1105/prod_lhxo
export BSG_REPRO_FULL_RUN=$BSG_REPRO_RUN BSG_BUILD_MERGED=
export BSG_INP=/resnick/groups/carnegie_poc/jingze/doface_1105/prep3d/bsgLI_lhx_exact_o.athinput
# IC overlays (Fig 3 right panel, r_photo,IC, Fe peak): L. Ma's profile (the run's he_ic_linhao_file)
# converted to the 7-column layout + his isothermal extension above rmax (tools/make_ic_linhao.py)
export BSG_IC=$W/ic/ic_linhao_7col.txt
export BSG_OPAC=/resnick/groups/carnegie_poc/jingze/bsgli_1014/files/rosseland_ma2026_x0.7_z0.008_padded.txt
export BSG_COST_PARTS=$W/parts.json BSG_RUN_LABEL=bsgLX-lhxo
export BSG_KE_LABEL="Linhao deck values + wrong v (vstar_nbe), GM 19.78, no FOFC, CFL 0.3, twin off, no leak mimic"
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
python3 $W/pipe/cmp_arms.py $PAGE/img/ours_arms.png > $W/pipe/cache/cmp_arms.txt 2>&1 || echo "cmp_arms FAILED"; tail -4 $W/pipe/cache/cmp_arms.txt
cp $W/pipe/index.html $PAGE/index.html
bash $W/pipe/update_page_repro.sh $PAGE
python3 $W/tools/dark_paper_fig4.py $PAGE   # 10-14: img/paper_fig4_dark.png (Fig 4 dark panel)
python3 $W/tools/inject_status.py $PAGE
