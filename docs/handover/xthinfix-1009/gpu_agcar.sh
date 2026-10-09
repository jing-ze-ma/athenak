#!/bin/bash
# xthinfix-1009 GPU part (H200): AG Car A production input (t = 0), 30 cycles, mode all vs beam,
# interleaved with repeats, on 2 H200 (2 ranks x 2 MeshBlocks, as the AG Car link scripts), then a
# 1-H200 point (4 MeshBlocks on one rank) for the scaling pair.  Adapt only the #SBATCH header and the
# srun/binding line to your link_agcarb.sh / link_agcarb_1g.sh (CUDA_VISIBLE_DEVICES=SLURM_LOCALID).
# usage: ATHENA_GPU=<he_star_m1 H200 binary of xthinfix-1009> IN_A=<files>/agcar_rcxA_ge_accel_st_mg_hr.athinput \
#        XTF_RUN=<run dir> REPO=<athenak checkout> BUNDLE=<this dir> sbatch --export=ALL gpu_agcar.sh
export OMP_NUM_THREADS=1
export HSA_XNACK=1; export HSA_NO_SCRATCH_RECLAIM=1   # job-env rule (no-op on CUDA)
: "${ATHENA_GPU:?}" "${IN_A:?}" "${XTF_RUN:?}" "${REPO:?}" "${BUNDLE:?}"; export REPO
S=$XTF_RUN/gpu.${SLURM_JOB_ID:-local}; mkdir -p $S
# the new keys must exist in the input file (a command-line key absent from a fresh input is FATAL)
IN=$S/agcar_A_xtf.athinput
python3 - "$IN_A" "$IN" <<'PY'
import sys
s = open(sys.argv[1]).read()
assert 'implicit_blend_xthin = 30' in s, 'not the hr production input'
s = s.replace('implicit_blend_xthin = 30', 'implicit_blend_xthin = 30\nimplicit_blend_xthin_mode = all\n'
              'implicit_blend_xthin_wmin = 0.0', 1)
open(sys.argv[2], 'w').write(s)
PY
O="time/nlim=30 time/ndiag=1 output5/dt=1e30"
run(){ D=$S/$1; NR=$2; shift 2; mkdir -p $D && cd $D || exit 1
  echo "== $(basename $D) start $(date +%T)"
  srun -n $NR --cpu-bind=cores $ATHENA_GPU -i $IN -d $D $O "$@" < /dev/null > run.log 2>&1
  echo "rc=$? fatal=$(grep -c FATAL run.log) $(grep 'cycle=' run.log | tail -n 1) $(date +%T)"; }
# 2 GPUs: interleaved arms
run all_1 2 rad_m1/implicit_blend_xthin_mode=all
run beam_1 2 rad_m1/implicit_blend_xthin_mode=beam
run all_2 2 rad_m1/implicit_blend_xthin_mode=all
run beam_2 2 rad_m1/implicit_blend_xthin_mode=beam
# 1 GPU scaling point (needs ~80 GB; same 4 MeshBlocks on one rank)
CUDA_VISIBLE_DEVICES=0 run beam_1g 1 rad_m1/implicit_blend_xthin_mode=beam
CUDA_VISIBLE_DEVICES=0 run all_1g 1 rad_m1/implicit_blend_xthin_mode=all
{
for a in all_1 beam_1 all_2 beam_2 beam_1g all_1g; do
  python3 - "$S/$a/run.log" "$a" <<'PY'
import re, sys
t = open(sys.argv[1]).read()
el = [(float(a), int(b)) for a, b in re.findall(r'elapsed=([0-9.e+-]+) cycle=([0-9]+)', t)]
pic = re.findall(r'Picard iterations mean=([0-9.e+-]+) max=([0-9.e+-]+) NON-CONVERGED=([0-9.e+-]+)', t)
inn = re.findall(r'inner iterations mean=([0-9.e+-]+)', t)
s = '?'
if len(el) > 6:
    a, b = el[5], el[-1]
    s = '%.3f s/cycle (cycles %d-%d)' % ((b[0] - a[0])/max(b[1] - a[1], 1), a[1], b[1])
print('%-8s %s  Picard mean/max/NONCONV %s  inner mean %s  FATAL %d' % (
    sys.argv[2], s, pic[-1] if pic else '?', inn[-1] if inn else '?', t.count('FATAL')))
PY
done
echo "-- bitwise repeats (hst+bin): all_1 vs all_2, beam_1 vs beam_2"
for p in "all_1 all_2" "beam_1 beam_2"; do set -- $p; n=0; d=0
  for f in $(cd $S/$1 && find . -name '*.hst' -o -name '*.bin' | sort); do n=$((n+1)); cmp -s $S/$1/$f $S/$2/$f || d=$((d+1)); done
  echo "$1 vs $2: $n files, $d differ"; done
echo "-- all vs beam by radius (max over angles of |a-b|/max|b|), last dump"
python3 $BUNDLE/ana/spdiff.py $S/all_1 $S/beam_1
} > $S/RESULTS_gpu.txt 2>&1
echo GPU_DONE $S/RESULTS_gpu.txt
