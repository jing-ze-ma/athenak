#!/bin/bash
# refresh.sh: rebuild the bsg_doface page (BSG Athena++-equivalent arm, doppler_vel = vstar_n; run doface_1105/prod_doface),
# CPU only, inside the held allocation 5719508 on hpc-sm-01-05 via srun --overlap (8 cores, 60 GB, nice; the 4 GPUs run the arm;
# no GPU use, no sbatch), then write the files to publish (index.html + media; no inlining, original resolution, USER 10-24).
# make_page.sh also runs pipe/cmp_arms.py (this arm vs the Linhao-exact production bsgli_1104/prod -> img/ours_arms.png).
# usage: setsid nohup bash /resnick/groups/carnegie_poc/jingze/doface_1105/page/refresh.sh [full] > /resnick/groups/carnegie_poc/jingze/doface_1105/page/refresh.out 2>&1 < /dev/null &
#        (full = redraw all frames; colour limits are fixed constants read off Ma+2026 since 11-06)
W=/resnick/groups/carnegie_poc/jingze/doface_1105/page; P=/resnick/home/jingze/ATHENAK/pages/bsg_doface
[ "$1" = full ] && export BSG_RECHOOSE=1 || export BSG_RECHOOSE=0
srun --overlap --jobid=5719508 -w hpc-sm-01-05 -N1 -n1 -c 4 --mem=60G nice -n 19 bash -l $W/make_page.sh $P > $W/make_page.out 2>&1
tail -3 $W/make_page.out
cd $P && { echo index.html; echo status.json; echo movie/n.json; ls img/*.png movie/*.mp4 movie/*/*.png 2>/dev/null; } | while read f; do [ -f "$f" ] && echo "$f"; done > $W/publish_files.txt
echo "publish $(wc -l < $W/publish_files.txt) files from $P (list: $W/publish_files.txt); total $(cd $P && du -cb $(cat $W/publish_files.txt) | tail -1 | cut -f1) B; largest: $(cd $P && du -b $(cat $W/publish_files.txt) | sort -n | tail -1)"
