#!/bin/bash
# refresh.sh: rebuild the bsg_lhx_none page (CPU only, inside the held allocation 5719508 on hpc-sm-01-05 via srun --overlap, 8 cores;
# no sbatch), then print the files to publish (index.html + media; no inlining, USER 10-24).
# usage: setsid nohup bash page_none/refresh.sh [full] > page_none/refresh.out 2>&1 &   (full = re-choose colour ranges, redraw all frames)
W=/resnick/groups/carnegie_poc/jingze/bsgli_1104/page_none; P=/resnick/home/jingze/ATHENAK/pages/bsg_lhx_none
[ "$1" = full ] && export BSG_RECHOOSE=1 || export BSG_RECHOOSE=0
bash /resnick/groups/carnegie_poc/jingze/bsgli_1104/run_cmp_mdot.sh; cp /resnick/groups/carnegie_poc/jingze/bsgli_1104/cmp_mdot.png $P/img/ours_mdot_cmp.png
srun --overlap --jobid=5719508 -w hpc-sm-01-05 -N1 -n1 -c 8 --mem=64G nice -n 10 bash -l $W/make_page.sh $P > $W/make_page.out 2>&1
tail -3 $W/make_page.out
cd $P && { echo index.html; echo status.json; echo movie/n.json; ls img/*.png movie/*.mp4 movie/*/*.png 2>/dev/null; } | while read f; do [ -f "$f" ] && echo "$f"; done > $W/publish_files.txt
echo "publish $(wc -l < $W/publish_files.txt) files from $P (list: $W/publish_files.txt); largest: $(cd $P && du -b $(cat $W/publish_files.txt) | sort -n | tail -1)"
