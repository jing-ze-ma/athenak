#!/bin/bash
# tools/ic_star regression + reproduction checks (login node is fine; ~20-40 min in total).
#   1. AG Car A deep70 (E1 production IC) from specs/agcar_A_deep70.yaml
#   2. hg145 ic2 (make_ic_agcar_mltpp.py --mode inward) from specs/hegiant_14p5M_hg145compat.yaml
# Output: $OUT/<name>/diag_<name>.txt (the 'regression:' block at the end).
set -e
cd "$(dirname "$0")/.."
PY=${PY:-/u/jinma/miniconda3/envs/magritte/bin/python3}
OUT=${OUT:-/raven/ptmp/jinma/ic_star_runs}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
$PY make_star_ic.py specs/agcar_A_deep70.yaml --out $OUT/agcar_A_deep70 \
    --compare /raven/ptmp/jinma/agcmodel_1012/ic/ic_agcar_A_ge_deep70.txt &
$PY make_star_ic.py specs/hegiant_14p5M_hg145compat.yaml --out $OUT/hegiant_14p5M_hg145compat \
    --compare /raven/ptmp/jinma/hegiant145_1015/ic2/ic_hg145_mltpp.txt &
wait
grep -A12 '^regression' $OUT/agcar_A_deep70/diag_agcar_A_deep70.txt
grep -A12 '^regression' $OUT/hegiant_14p5M_hg145compat/diag_hegiant_14p5M_hg145compat.txt
