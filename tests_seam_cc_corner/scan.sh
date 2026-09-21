#!/bin/bash
# the static corner/edge/face seam-halo scan, all layouts.  $1 = binary.
# mesh/nx1=16 meshblock/nx1=8 gives the TWO radial blocks the x1x2x3 corner slots need.
B=${1:-./athena_after}
I=../inputs/tests/cubed_sphere_raddiff.athinput
for n in 32 64; do
  for mb in $n $((n/2)) $((n/4)); do
    echo "=== n=$n mb=$mb"
    $B -i $I -d out mesh/nx1=16 meshblock/nx1=8 mesh/nx2=$n mesh/nx3=$n \
       meshblock/nx2=$mb meshblock/nx3=$mb time/nlim=0 problem/seam_halo_scan=1 \
       output1/dt=1e9 2>&1 | grep -E "GLOBAL MAX|max ="
    rm -rf out
  done
done
