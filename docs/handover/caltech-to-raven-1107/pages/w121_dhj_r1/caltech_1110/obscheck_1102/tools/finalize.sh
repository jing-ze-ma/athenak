#!/bin/bash
# waits for a free agent GPU, runs all em3d passes, bands, plots, takeaways, page rebuild (no publish)
O=/resnick/groups/carnegie_poc/jingze/obscheck_1102; cd $O; exec > $O/finalize.log 2>&1
until tools/anygpu.sh tools/run_em.sh; do echo "$(date) no GPU yet"; sleep 120; done
srun --jobid=5274196 --overlap -n1 -c2 bash -c "venv/bin/python tools/bands.py w121 out/bands_w121.npz out/em_w121_0059[6-9].npz out/em_w121_00600.npz && venv/bin/python tools/bands.py w76 out/bands_w76.npz out/em_w76_0059[6-9].npz out/em_w76_00600.npz && venv/bin/python tools/obs_plots.py && venv/bin/python tools/takeaways.py"
/resnick/home/jingze/ATHENAK/pages/dhj_r1/refresh.sh
echo FINALIZE DONE $(date)
