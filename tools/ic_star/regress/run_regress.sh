#!/bin/bash
# tools/ic_star regressions (login node is fine but slow; ~30-60 min in total, runs in parallel):
#   1. mode global, agcar_legacy: AG Car A deep70 (E1 production IC) -- max rel diff (round-off)
#   2. mode global, scheme mltpp_inward: hg145 ic2 (make_ic_agcar_mltpp.py) -- bit for bit
#   3. mode mesa: 14p5M vs MESA (m_env, R_ph, CZs) and vs hg145 ic3 (same method, own builder)
#   4. mode mesa: BSG vs MESA (R_ph shot, m_env)
set -e
cd "$(dirname "$0")/.."
PY=${PY:-/u/jinma/miniconda3/envs/magritte/bin/python3}
OUT=${OUT:-/raven/ptmp/jinma/ic_star_runs}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
$PY make_star_ic.py specs/agcar_A_deep70.yaml --no-scan --out $OUT/agcar_A_deep70 \
    --compare /raven/ptmp/jinma/agcmodel_1012/ic/ic_agcar_A_ge_deep70.txt > $OUT/regress_agcar.log 2>&1 &
$PY make_star_ic.py specs/hegiant_14p5M_hg145compat.yaml --no-scan --out $OUT/hegiant_14p5M_hg145compat \
    --compare /raven/ptmp/jinma/hegiant145_1015/ic2/ic_hg145_mltpp.txt > $OUT/regress_hg145c.log 2>&1 &
$PY make_star_ic.py specs/hegiant_14p5M.yaml --out $OUT/hegiant_14p5M > $OUT/regress_145.log 2>&1 &
$PY make_star_ic.py specs/bsg.yaml --out $OUT/bsg > $OUT/regress_bsg.log 2>&1 &
wait
grep -A12 '^regression' $OUT/agcar_A_deep70/diag_agcar_A_deep70.txt
grep -A12 '^regression' $OUT/hegiant_14p5M_hg145compat/diag_hegiant_14p5M_hg145compat.txt
for n in hegiant_14p5M bsg; do
  echo "== mode mesa $n"; grep -E 'R_ph ours|m\(x1min\) ours|F_rad\(our kappa\)|^    [0-9.]+-[0-9.]+ Rsun' $OUT/$n/diag_$n.txt
done
$PY regress/cmp_radii.py $OUT/hegiant_14p5M/ic_hegiant_14p5M.txt \
    /raven/ptmp/jinma/hegiant145_1015/ic3/ic_hg145_mesaP.txt 6.5 8 10 19.4 30 50 80 100 108
