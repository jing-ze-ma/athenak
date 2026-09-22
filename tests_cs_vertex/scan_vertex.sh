#!/bin/bash
# GATE (b): the STATIC cube-vertex ghost error, cell-centred and face-centred.
#
#   ./scan_vertex.sh <binary>
#
# CELL-CENTRED: cs_test iprob 15 (a spherical harmonic of the direction cosines, whose
# exact value in a ghost cell is known from the ghost's OWN chart-continued (xi,eta)),
# run with time/nlim = 0 so exactly one boundary exchange is measured.  The scan's
# category bins already separate "edge, CUBE VERTEX" (the tangential vertex block) and
# "corner, CUBE VRTX" (the same block in a RADIAL ghost layer, which needs two radial
# MeshBlocks -- hence mesh/nx1=16 meshblock/nx1=8).  Arms: mesh/cs_vertex_fill_cc.
#
# FACE-CENTRED: cs_test iprob 11 with seam_halo_scan_fc, arms mesh/cs_vertex_fill --
# the switch that already exists -- to confirm the recorded extrapolate-vs-sample gain.
set -u
B=${1:-./athena_new}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE" || exit 1
mkdir -p out

echo "######## CELL-CENTRED (iprob 15), arms: cs_vertex_fill_cc = false | true"
for n in 32 64; do
  for mb in $n $((n/2)); do
    for v in false true; do
      echo "=== cc n=$n mb=$mb cs_vertex_fill_cc=$v"
      $B -i in_scan_cc.athinput -d out mesh/nx1=16 meshblock/nx1=8 \
         mesh/nx2=$n mesh/nx3=$n meshblock/nx2=$mb meshblock/nx3=$mb \
         time/nlim=0 problem/seam_halo_scan=1 mesh/cs_vertex_fill_cc=$v \
         2>&1 | grep -E "GLOBAL MAX|max ="
      rm -rf out; mkdir -p out
    done
  done
done

echo "######## FACE-CENTRED (iprob 11), arms: cs_vertex_fill = false | true"
for n in 32 64; do
  for mb in $n $((n/2)); do
    for v in false true; do
      echo "=== fc n=$n mb=$mb cs_vertex_fill=$v"
      $B -i in_scan_fc.athinput -d out mesh/nx1=16 meshblock/nx1=8 \
         mesh/nx2=$n mesh/nx3=$n meshblock/nx2=$mb meshblock/nx3=$mb \
         time/nlim=0 problem/seam_halo_scan_fc=1 mesh/cs_vertex_fill=$v \
         2>&1 | grep -E "GLOBAL MAX|max ="
      rm -rf out; mkdir -p out
    done
  done
done
rm -rf out
