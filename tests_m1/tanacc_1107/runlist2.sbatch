#!/bin/bash -l
#SBATCH -p apudev
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --constraint=apu
#SBATCH --cpus-per-task=24
#SBATCH --time 00:15:00
#SBATCH -o /viper/ptmp2/jinma/openitems_1107/tanacc/logs/%x.%j.out
# usage: sbatch -J name runlist.sbatch <listfile> <runs_root> [extra keys...]
# runs each deck name of <listfile> (inp/<name>.athinput) in <runs_root>/<name>, two at a time (1 GPU each)
module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0
export OMP_NUM_THREADS=1 HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1
T=/viper/ptmp2/jinma/openitems_1107/tanacc
BIN=$T/bin/athena_tanfix; MD5=$(cat $T/bin/athena_tanfix.md5)
[ "$(md5sum $BIN | cut -c1-32)" = "$MD5" ] || { echo "md5 mismatch"; exit 1; }
LIST=$1; ROOT=$2; shift 2; EXTRA="$*"
go() {
  local n=$1 d=$ROOT/$1; mkdir -p $d; cd $d || return 1
  cp -f $T/inp2/$n.athinput $d/
  echo "$BIN -i $T/inp2/$n.athinput $EXTRA -d $d" > cmd.$SLURM_JOB_ID.txt
  local t0=$(date +%s)
  srun -n 1 --exact --gpus-per-task=1 --cpus-per-task=24 $BIN -i $T/inp2/$n.athinput $EXTRA -d $d \
    < /dev/null > run.log 2>&1; echo "rc=$? wall=$(( $(date +%s) - t0 ))" >> run.log; }
echo "job $SLURM_JOB_ID $SLURM_JOB_NODELIST $(date +%F_%T) list $LIST extra '$EXTRA'"
mapfile -t N < $LIST
for ((i=0; i<${#N[@]}; i+=2)); do
  go ${N[i]} & [ -n "${N[i+1]}" ] && go ${N[i+1]} &
  wait
  echo "$(date +%T) ${N[i]} $(tail -1 $ROOT/${N[i]}/run.log) | ${N[i+1]} $(tail -1 $ROOT/${N[i+1]}/run.log 2>/dev/null)"
done
echo "end $(date +%T)"
