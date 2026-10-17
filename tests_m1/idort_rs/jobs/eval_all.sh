#!/bin/bash
# idort-rs-1017: evaluate every finished rsx run (rc.txt = 0) -> eval/rsx.txt, rsx.json, rsx_table.md
W=${W:-/raven/ptmp/jinma/idortrs_1017}; H=$(cd $(dirname $0)/.. && pwd); P=${PY:-/u/jinma/miniconda3/envs/magritte/bin/python}
export OPENBLAS_NUM_THREADS=1
ds=$(for f in $(find $W/runs -name rc.txt | sort); do [ "$(cat $f)" = 0 ] && dirname $f; done)
RS_JSON=$W/eval/rsx.json $P $H/rs_eval.py $ds > $W/eval/rsx.txt 2>&1
$P $H/mk_table.py $W/eval/rsx.json > $W/eval/rsx_table.md
echo "evaluated $(echo $ds | wc -w) runs"
