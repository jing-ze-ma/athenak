#!/bin/bash
# tanfix-1109 cost on one Raven node (4 A100): GPU q runs case q's arms interleaved off on on_tw on_hd, x2
module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0
export OMP_NUM_THREADS=1 HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1
R=/raven/ptmp/jinma/openitems_1107/tanacc
BIN=$R/bin/athena_tanfix; MD5=$(cat $R/bin/athena_tanfix.md5)
[ "$(md5sum $BIN | cut -c1-32)" = "$MD5" ] || { echo "md5 mismatch"; exit 1; }
sq() { local c=$1
  for rep in 1 2; do for a in off on on_tw on_hd; do
    n=${c}_$a; d=$R/timing/$n.r$rep; mkdir -p $d; cd $d || return
    python3 /raven/ptmp/jinma/pool_1107/preflight.py --md5 $MD5 $BIN -i $R/inp2/$n.athinput > /dev/null || { echo "$n PREFLIGHT"; continue; }
    t0=$(date +%s.%N)
    srun -n 1 --exact --gpus-per-task=1 --cpus-per-task=18 --mem=100G $BIN -i $R/inp2/$n.athinput -d $d < /dev/null > run.log 2>&1
    echo "$n rep$rep rc=$? wall=$(echo "$(date +%s.%N) - $t0" | bc)"
  done; done; }
mkdir -p $R/timing
echo "job $SLURM_JOB_ID $SLURM_JOB_NODELIST"
for c in sp_hi_hot cs_hi_hot sp_hi_shd cs_hi_shd; do sq $c > $R/timing/$c.$SLURM_JOB_ID.txt & done; wait
cat $R/timing/*.$SLURM_JOB_ID.txt
