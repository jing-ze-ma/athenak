#!/bin/bash
# Rebuild the Linhao BSG page (CPU debug job, waits for it), then print the JSON list of files to republish
# (fixed set + new movie frames). Usage: bash upd_page.sh   (run in background; then publish index.html with root=$P and that list)
W=/resnick/groups/carnegie_poc/jingze/bsgli_1014/pagework_li; P=/resnick/home/jingze/ATHENAK/pages/bsg_li_caltech
bash $W/merge_run.sh; bash $W/make_page.sh $P > $W/make_page.login.out 2>&1
J=$W/make_page.login.out; tail -3 $J
bash /resnick/groups/carnegie_poc/jingze/vis3d/refresh.sh bsg_li   # 3D movie: submits a CPU render job only if new dumps (non-blocking)
cd $P && { echo status.json; echo movie/n.json; ls img/*.png; [ -f new_files.txt ] && grep -E '^movie/' new_files.txt; ls movie/*.mp4 2>/dev/null; ls movie3d/*.mp4 movie3d/*.png 2>/dev/null; } \
  | sort -u | while read f; do [ -f "$f" ] && echo "$f"; done | python3 -c 'import sys,json; print(json.dumps([l.strip() for l in sys.stdin]))'
