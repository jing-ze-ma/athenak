#!/bin/bash
# Re-verification of the FC cube-vertex claim in tests_cs_vertex/README.md.
# Runs the static seam-halo scan (iprob 11, seam_halo_scan_fc) at n = 32/64/128,
# 1x1 and 2x2 MeshBlocks/panel, both cs_vertex_fill = false/true, and BOTH
# faces_from_potential = true (as the original scan_vertex.sh ran it, via the
# input file's own default) and = false (as the function's own docstring says
# it must be run, to remove the O(h^2) discretisation offset of the active
# faces themselves from the comparison).
set -u
B=${1:-./athena_new}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE" || exit 1
mkdir -p out

for fp in true false; do
  echo "########################################## faces_from_potential=$fp"
  for n in 32 64 128; do
    for mb in $n $((n/2)); do
      for v in false true; do
        echo "=== fc n=$n mb=$mb cs_vertex_fill=$v faces_from_potential=$fp"
        $B -i in_scan_fc.athinput -d out mesh/nx1=8 meshblock/nx1=8 \
           mesh/nx2=$n mesh/nx3=$n meshblock/nx2=$mb meshblock/nx3=$mb \
           time/nlim=0 problem/seam_halo_scan_fc=1 mesh/cs_vertex_fill=$v \
           problem/faces_from_potential=$fp \
           2>&1 | grep -E "GLOBAL MAX|CUBE VERTEX"
        rm -rf out; mkdir -p out
      done
    done
  done
done
rm -rf out
