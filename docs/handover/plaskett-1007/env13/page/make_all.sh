#!/bin/bash
# One command: figures + numbers (mkfigs13.py) for q13_s1 and the q12 reference arms, then the
# page (build_page13.py) -> index.html (+ figs/) and index_embedded.html (single file).
#
#   ./make_all.sh                       # production: run/q13_s1 -> page13/index.html, page13/figs/
#   ./make_all.sh --standin <run dir>   # pipeline test: <run dir> plays q13_s1 -> page13/standin_test/
#   extra args after those (e.g. --state running, --force) are passed on: --state to the
#   builder, --force to mkfigs (redo maps/edge caches)
set -euo pipefail
P=/work/nvme/bivj/jma20/plaskett_1007/test1007/page13
T=/work/nvme/bivj/jma20/plaskett_1007/test1007
R12=/work/nvme/bivj/jma20/plaskett_1007/run
PY=/u/jma20/ATHENAK/venv/bin/python          # numpy, scipy, matplotlib
RUN=$T/run/q13_s1; OUT=$P; SI=""
STATE=auto; FORCE=""
while [ $# -gt 0 ]; do
  case $1 in
    --standin) RUN=$(readlink -f "$2"); OUT=$P/standin_test; SI=--standin; shift 2;;
    --state) STATE=$2; shift 2;;
    --force) FORCE=--force; shift;;
    *) echo "unknown arg $1"; exit 2;;
  esac
done
mkdir -p "$OUT/figs"
echo "== mkfigs13: q13_s1 <- $RUN ${SI:+(STAND-IN)}"
$PY $P/mkfigs13.py --run "$RUN" $R12/q12_s1 $R12/q12f_s1 --names q13_s1 q12_s1 q12f_s1 \
    --out "$OUT/figs" --vis $P/vis $FORCE \
    --envonly t2_sphere=$T/t2_sphere t2_eq=$T/t2_eq t4_env13=$T/t4_env13
echo "== build_page13"
$PY $P/build_page13.py --run "$RUN" --figs "$OUT/figs" --out "$OUT" --state "$STATE" $SI
$PY $P/build_page13.py --run "$RUN" --figs "$OUT/figs" --out "$OUT" --state "$STATE" $SI --embed
echo "== done: $OUT/index.html (+ figs/*.png) and $OUT/index_embedded.html"
ls -la "$OUT"/index*.html "$OUT"/figs/*.png | awk '{print $5, $9}'
