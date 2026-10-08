#!/bin/bash
# usage: bash run_batch1.sh <binary>   (from this directory; submits 5 jobs, ~4.5 GPU-h at most)
BIN=$1; H=$(cd "$(dirname "$0")" && pwd)
mkdir -p /work/nvme/bivj/jma20/delta_1008/spb2
P=$(sbatch --parsable $H/spb2_delta.sh $BIN A prep $H/armsA0.txt 236 00:13:00)
echo "prepA $P"
for a in A1 A2 A3; do
  echo "$a $(sbatch --parsable --dependency=afterok:$P $H/spb2_delta.sh $BIN A arms $H/arms$a.txt 276 00:02:40)"
done
echo "B1 $(sbatch --parsable $H/spb2_delta.sh $BIN B arms $H/armsB1.txt 40 00:03:00)"
