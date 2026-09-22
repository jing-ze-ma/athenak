#!/bin/bash
# Knob sweep on the cs field loop at n = 32, t = 2.0.  Two concurrent solvers.
set -u
H="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$H/out" || exit 1
NPROC=${NPROC:-2}
run() {
  local tag=$1; shift
  [ -s "$H/out/${tag}.user.hst" ] && { echo "skip $tag"; return 0; }
  "$H/athena_diag" -i "$H/loop_hist.athinput" -d . job/basename="$tag" \
    mesh/nx1=32 mesh/nx2=32 mesh/nx3=32 \
    meshblock/nx1=32 meshblock/nx2=32 meshblock/nx3=32 \
    problem/conv_nband=4 "$@" > "$H/logs/${tag}.log" 2>&1
}
while read -r tag rest; do
  [ -z "$tag" ] && continue
  while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 5; done
  echo "run $tag $rest"
  run "$tag" $rest &
done <<'ARMS'
k_bsemf mhd/bs_emf=true
k_vfillcc mesh/cs_vertex_fill_cc=true
k_ppmx mhd/reconstruct=ppmx
k_fofc mhd/fofc=true
k_nolowbeta mhd/cs_lowbeta_fallback=0.0
k_wbsrc mhd/cs_wellbalanced_src=true
k_novfill mesh/cs_vertex_fill=false
k_wenoz mhd/reconstruct=wenoz
ARMS
wait
