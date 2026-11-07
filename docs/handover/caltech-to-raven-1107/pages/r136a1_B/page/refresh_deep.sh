#!/bin/bash
# refresh.sh: rebuild the R136a1 page from NEW dumps only (CPU, inside the held allocation 5544871 on hpc-sm-02-03 via
# srun --overlap, 8 cores, nice 10; no GPU, no sbatch), then write publish_files_deep.txt (index.html + media; no inlining).
# usage: setsid nohup bash /resnick/groups/carnegie_poc/jingze/r136a1_1104/page/refresh_deep.sh [full] > .../page/refresh_deep.out 2>&1 < /dev/null &
W=/resnick/groups/carnegie_poc/jingze/r136a1_1104/page; P=/resnick/home/jingze/ATHENAK/pages/r136a1_deep
JID=${JID:-5544871}; ND=${ND:-hpc-sm-02-03}
[ "$1" = full ] && F=--full || F=
srun --overlap --jobid=$JID -w $ND -N1 -n1 -c 4 --mem=16G nice -n 19 bash -l $W/make_page_deep.sh $P $F > $W/make_page_deep.out 2>&1
tail -2 $W/make_page_deep.out
cd $P && { echo index.html; echo status.json; echo movie/n.json; ls img/*.png movie/eq/*.png movie/me/*.png 2>/dev/null; } \
  | while read f; do [ -f "$f" ] && echo "$f"; done > $W/publish_files_deep.txt
cp $W/publish_files_deep.txt $P/publish_files.txt
echo "publish $(wc -l < $W/publish_files_deep.txt) files from $P (list: $W/publish_files_deep.txt); largest: $(cd $P && du -b $(cat $W/publish_files_deep.txt) | sort -n | tail -1)"
