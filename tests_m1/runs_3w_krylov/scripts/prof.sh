#!/bin/bash -l
# rocprofv3 kernel + HIP runtime trace; launches and host syncs per cycle from the
# difference of an nlim=20 and an nlim=10 run.  usage: sbatch prof.sh <tag> <armfile>
#SBATCH -o /viper/ptmp2/jinma/krylov_0923/runs/prof.out.%j
#SBATCH -e /viper/ptmp2/jinma/krylov_0923/runs/prof.err.%j
#SBATCH -J m1kprof
#SBATCH -p apudev
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/krylov_0923
EXE=$W/bin/athena_new_gpu
TAG=$1; ARMF=$2
md5sum $EXE
mapfile -t ARMS < <(grep -v '^#' $W/scripts/$ARMF | grep -v '^ *$')
for line in "${ARMS[@]}"; do
  read -r name np inp ov <<< "$line"
  for nl in 10 20; do
    d=$W/runs/$TAG/${name}_n$nl; rm -rf $d; mkdir -p $d; cd $d
    srun -n $np --ntasks-per-node=2 bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec rocprofv3 --kernel-trace --hip-runtime-trace --output-format csv -d $PWD/prof -o rank_$SLURM_PROCID -- "$0" "$@"' \
      $EXE -i $W/gpu/$inp -d $d -t 00:03:00 time/nlim=$nl $ov > run.log 2> run.err
    echo "#### ${name}_n$nl rc=$? $(date +%T)"
  done
done
python3 $W/scripts/profsum.py $W/runs/$TAG
