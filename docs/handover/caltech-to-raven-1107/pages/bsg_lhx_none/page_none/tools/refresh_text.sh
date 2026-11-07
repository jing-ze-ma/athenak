#!/bin/bash
# light re-run of the text steps only (gen_status -> compare57 -> inject), login node OK
W=/resnick/groups/carnegie_poc/jingze/bsgli_1014/pagework_li; PAGE=${1:-/resnick/home/jingze/ATHENAK/pages/bsg_li_caltech}
eval "$(sed -n '/^source \/etc\/profile/,/^export MPLCONFIGDIR/p' $W/make_page.sh)"
cd $W/pipe
nice python3 gen_status.py $PAGE >/dev/null || exit 1
nice python3 compare57.py $PAGE > cache/cmp57_agg.log 2>&1 && cp out/ours_cmp57.png $PAGE/img/
python3 $W/tools/inject_status.py $PAGE
