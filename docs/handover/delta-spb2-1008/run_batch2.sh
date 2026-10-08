#!/bin/bash
# usage: bash run_batch2.sh <binary of batch 2>   (needs $W/prepA/rst from batch 1; 3 jobs, ~2.7 GPU-h)
BIN=$1; H=$(cd "$(dirname "$0")" && pwd)
ls /work/nvme/bivj/jma20/delta_1008/spb2/prepA/rst/*.rst || { echo "no prepA restart (batch 1)"; exit 1; }
for a in A4 A5; do
  echo "$a $(sbatch --parsable $H/spb2_delta.sh $BIN A arms $H/arms$a.txt 276 00:02:40)"
done
echo "B2 $(sbatch --parsable $H/spb2_delta.sh $BIN B arms $H/armsB2.txt 40 00:03:00)"
