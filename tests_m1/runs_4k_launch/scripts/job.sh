#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/launch_0923/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/launch_0923/runs/log.err.%j
#SBATCH -J m1launch
#SBATCH --constraint="apu"
#SBATCH --cpus-per-task=24
# usage: sbatch -p apudev --gres=gpu:N --ntasks-per-node=N [-N nodes] -t T job.sh <tag> <armfile>
#   env: BINS="name=path ..." (arm line picks one by its 3rd field), NREP (default 3),
#        NLIM (default 120), PROF=1 (rocprofv3 kernel + HIP runtime trace, one repeat)
# arm file lines: <name> <np> <binname> <input> [overrides...]
# repeats: odd in file order, even reversed (interleaved arms).
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/launch_0923
TAG=$1; ARMF=$2; [ -f $W/scripts/env_$TAG.sh ] && source $W/scripts/env_$TAG.sh
NREP=${NREP:-3}; NLIM=${NLIM:-120}; PROF=${PROF:-0}
declare -A B
for kv in $BINS; do B[${kv%%=*}]=${kv#*=}; md5sum ${kv#*=}; done
mapfile -t AL < <(grep -v '^#' $W/scripts/$ARMF | grep -v '^ *$')
[ "$PROF" = 1 ] && NREP=1
for r in $(seq 1 $NREP); do
  if [ $((r % 2)) -eq 0 ]; then ORD=$(seq $((${#AL[@]}-1)) -1 0); else ORD=$(seq 0 $((${#AL[@]}-1))); fi
  for x in $ORD; do
    read -r name np bn inp ov <<< "${AL[$x]}"
    d=$W/runs/$TAG/${name}_$r; rm -rf $d; mkdir -p $d; cd $d
    echo "#### ${name}_$r $(date +%T)"
    if [ "$PROF" = 1 ]; then
      srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec rocprofv3 --kernel-trace --hip-runtime-trace --output-format csv -d $PWD/prof -o rank_$SLURM_PROCID -- "$0" "$@"' \
        ${B[$bn]} -i $W/inp/$inp.athinput -d . -t 00:03:00 time/nlim=$NLIM $ov > run.log 2> run.err
    else
      srun -n $np bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
        ${B[$bn]} -i $W/inp/$inp.athinput -d . -t 00:03:00 time/nlim=$NLIM $ov > run.log 2> run.err
    fi
    echo "rc=$? $(date +%T)"
  done
done
