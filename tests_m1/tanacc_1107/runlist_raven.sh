#!/bin/bash
# tanfix-1109 Raven runner (pool task or gpudev): bash runlist_raven.sh <listfile> <runs_root> [extra keys]
# runs the decks inp2/<name>.athinput 4 at a time, 1 A100 each; preflight per run; md5 check
module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0
export OMP_NUM_THREADS=1 HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1
R=/raven/ptmp/jinma/openitems_1107/tanacc
BIN=$R/bin/athena_tanfix; MD5=$(cat $R/bin/athena_tanfix.md5)
[ "$(md5sum $BIN | cut -c1-32)" = "$MD5" ] || { echo "md5 mismatch"; exit 1; }
LIST=$1; ROOT=$2; shift 2; EXTRA="$*"
go() {
  local n=$1 d=$ROOT/$1; mkdir -p $d; cd $d || return 1
  python3 /raven/ptmp/jinma/pool_1107/preflight.py --md5 $MD5 $BIN -i $R/inp2/$n.athinput $EXTRA > preflight.txt 2>&1 \
    || { echo "rc=PREFLIGHT" > run.log; return 1; }
  cp -f $R/inp2/$n.athinput $d/
  echo "$BIN -i $R/inp2/$n.athinput $EXTRA -d $d" > cmd.${SLURM_JOB_ID}.txt
  local t0=$(date +%s)
  srun -n 1 --exact --gpus-per-task=1 --cpus-per-task=18 --mem=100G $BIN -i $R/inp2/$n.athinput $EXTRA -d $d \
    < /dev/null > run.log 2>&1; echo "rc=$? wall=$(( $(date +%s) - t0 ))" >> run.log; }
echo "job $SLURM_JOB_ID $SLURM_JOB_NODELIST $(date +%F_%T) list $LIST extra '$EXTRA'"
mapfile -t N < $LIST
for ((i=0; i<${#N[@]}; i+=4)); do
  for q in 0 1 2 3; do [ -n "${N[i+q]}" ] && go ${N[i+q]} & done
  wait
  for q in 0 1 2 3; do [ -n "${N[i+q]}" ] && echo "$(date +%T) ${N[i+q]} $(tail -1 $ROOT/${N[i+q]}/run.log)"; done
done
echo "end $(date +%T)"
