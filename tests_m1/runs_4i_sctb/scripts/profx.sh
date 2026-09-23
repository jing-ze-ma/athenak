#!/bin/bash -l
# rocprofv3 kernel trace of the SC sweep.  usage: sbatch prof.sh <tag> <armfile>
#SBATCH -o /viper/ptmp2/jinma/sctb_0923/gpu/runs/prof.out.%j
#SBATCH -e /viper/ptmp2/jinma/sctb_0923/gpu/runs/prof.err.%j
#SBATCH -J sctb_prof
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/sctb_0923/gpu
EXE=$W/bin/${3:-athena_new_gpu}
TAG=$1; ARMF=$2
md5sum $EXE
mapfile -t ARMS < <(grep -v '^#' $W/$ARMF | grep -v '^ *$')
for line in "${ARMS[@]}"; do
  read -r name np inp ov <<< "$line"
  d=$W/runs/$TAG/$name; rm -rf $d; mkdir -p $d; cd $d
  srun -n $np --ntasks-per-node=2 bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec rocprofv3 --kernel-trace --output-format csv -d $PWD/prof -o rank_$SLURM_PROCID -- "$0" "$@"' \
    $EXE -i $W/$inp -d $d -t 00:03:00 time/nlim=12 $ov > run.log 2> run.err
  echo "#### $name rc=$? $(date +%T)"; grep "SC seconds" run.log
done
