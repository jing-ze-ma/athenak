#!/bin/bash
# xthinfix-1009 STEEP GPU point (H200): AG Car A hr production input, t = 0, 30 cycles, 2 H200
# (2 ranks x 2 MeshBlocks as your AG Car link scripts), arms interleaved:
#   base  = rt-integration 5b304cf9 binary, the production input unchanged
#   all   = xthinfix-1009 21795e2a binary, mode all (must be bin-data identical to base)
#   steep = mode steep (X0 30), steepp = mode steep + implicit_hr_recon = plm
# Adapt only the #SBATCH header and the srun / CUDA_VISIBLE_DEVICES line to your link scripts.
# usage: ATHENA_BASE=<5b304cf9 he H200> ATHENA_GPU=<21795e2a he H200> IN_A=<files>/agcar_rcxA_ge_accel_st_mg_hr.athinput \
#        XTF_RUN=<run dir> REPO=<xthinfix-1009 checkout> BUNDLE=<docs/handover/xthinfix-1009> sbatch --export=ALL gpu_steep.sh
export OMP_NUM_THREADS=1
export HSA_XNACK=1; export HSA_NO_SCRATCH_RECLAIM=1   # job-env rule (no-op on CUDA)
: "${ATHENA_BASE:?}" "${ATHENA_GPU:?}" "${IN_A:?}" "${XTF_RUN:?}" "${REPO:?}" "${BUNDLE:?}"; export REPO
S=$XTF_RUN/gpu_steep.${SLURM_JOB_ID:-local}; mkdir -p $S
IN=$S/agcar_A_xtf.athinput
python3 - "$IN_A" "$IN" <<'PY'
import sys
s = open(sys.argv[1]).read()
assert 'implicit_blend_xthin = 30' in s, 'not the hr production input'
s = s.replace('implicit_blend_xthin = 30', 'implicit_blend_xthin = 30\nimplicit_blend_xthin_mode = all\n'
              'implicit_blend_xthin_wmin = 0.0\nimplicit_hr_recon = dc\nimplicit_hr_recon_fresh = 2\n'
              'implicit_hr_damp = 0.0', 1)
open(sys.argv[2], 'w').write(s)
PY
O="time/nlim=30 time/ndiag=1 output5/dt=1e30"
run(){ D=$S/$1; X=$2; I=$3; shift 3; mkdir -p $D && cd $D || exit 1
  echo "== $(basename $D) start $(date +%T)"
  srun -n 2 --cpu-bind=cores $X -i $I -d $D $O "$@" < /dev/null > run.log 2>&1
  echo "rc=$? fatal=$(grep -c FATAL run.log) $(grep 'cycle=' run.log | tail -n 1) $(date +%T)"; }
run base_1 $ATHENA_BASE $IN_A
run all_1 $ATHENA_GPU $IN rad_m1/implicit_blend_xthin_mode=all
run steep_1 $ATHENA_GPU $IN rad_m1/implicit_blend_xthin_mode=steep
run steepp_1 $ATHENA_GPU $IN rad_m1/implicit_blend_xthin_mode=steep rad_m1/implicit_hr_recon=plm
run base_2 $ATHENA_BASE $IN_A
run steep_2 $ATHENA_GPU $IN rad_m1/implicit_blend_xthin_mode=steep
{
for a in base_1 all_1 steep_1 steepp_1 base_2 steep_2; do
  python3 - "$S/$a/run.log" "$a" <<'PY'
import re, sys
t = open(sys.argv[1]).read()
el = [(float(a), int(b)) for a, b in re.findall(r'elapsed=([0-9.e+-]+) cycle=([0-9]+)', t)]
pic = re.findall(r'Picard iterations mean=([0-9.e+-]+) max=([0-9.e+-]+) NON-CONVERGED=([0-9.e+-]+)', t)
inn = re.findall(r'inner iterations mean=([0-9.e+-]+)', t)
pos = re.findall(r'plm positivity fallbacks=([0-9.e+-]+)', t)
s = '?'
if len(el) > 6:
    a, b = el[5], el[-1]
    s = '%.3f s/cycle (cycles %d-%d)' % ((b[0] - a[0])/max(b[1] - a[1], 1), a[1], b[1])
print('%-9s %s  Picard mean/max/NONCONV %s  inner mean %s  plm pos-fallbacks %s  FATAL %d' % (
    sys.argv[2], s, pic[-1] if pic else '?', inn[-1] if inn else '?', pos[-1] if pos else '-',
    t.count('FATAL')))
PY
done
echo "-- bin data base_1 vs all_1 (after the parameter header <par_end>), and repeats"
for p in "base_1 all_1" "base_1 base_2" "steep_1 steep_2"; do set -- $p; n=0; d=0
  for f in $(cd $S/$1 && find . -name '*.bin' | sort); do n=$((n+1))
    oa=$(grep -abo '<par_end>' $S/$1/$f | head -1 | cut -d: -f1); ob=$(grep -abo '<par_end>' $S/$2/$f | head -1 | cut -d: -f1)
    cmp -s -i $oa:$ob $S/$1/$f $S/$2/$f || d=$((d+1)); done
  echo "$1 vs $2: $n bin files, $d differ in data"; done
echo "-- all vs steep / steep+plm by radius (max over angles |a-b|/max|b|), last dump"
python3 $BUNDLE/ana/spdiff.py $S/all_1 $S/steep_1
python3 $BUNDLE/ana/spdiff.py $S/all_1 $S/steepp_1
} > $S/RESULTS_gpu_steep.txt 2>&1
echo GPU_DONE $S/RESULTS_gpu_steep.txt
