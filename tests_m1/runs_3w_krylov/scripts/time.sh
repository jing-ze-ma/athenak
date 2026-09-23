#!/bin/bash -l
# M1 cycle cost on the GPU: ONE binary, interleaved arms, repeat a in order and repeat b
# in reverse order.  usage: sbatch [-p ..] [-N n] time.sh <tag> <armfile>
# armfile lines: <name> <nranks> <input> [overrides...]
#SBATCH -o /viper/ptmp2/jinma/krylov_0923/runs/time.out.%j
#SBATCH -e /viper/ptmp2/jinma/krylov_0923/runs/time.err.%j
#SBATCH -J m1krylov
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/krylov_0923
EXE=$W/bin/athena_new_gpu
TAG=$1; ARMF=$2
md5sum $EXE $W/bin/athena_ref_gpu
mapfile -t ARMS < <(grep -v '^#' $W/scripts/$ARMF | grep -v '^ *$')
run() {
  local line="$1" rep=$2
  read -r name np inp ov <<< "$line"
  local d=$W/runs/$TAG/${name}_$rep; rm -rf $d; mkdir -p $d; cd $d
  echo "#### ${name}_$rep np=$np inp=$inp ov=$ov node=$(hostname) $(date +%T)"
  local exe=$EXE; [[ $name == head_* ]] && exe=$W/bin/athena_ref_gpu
  srun -n $np --ntasks-per-node=2 bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
    $exe -i $W/gpu/$inp -d $d -t 00:02:30 time/nlim=100 $ov > run.log 2> run.err
  echo "rc=$? $(date +%T)"; grep -E "NON-CONVERGED|inner iterations|breakdowns|halo_mpi" run.log
}
for a in "${ARMS[@]}"; do run "$a" a; done
for ((i=${#ARMS[@]}-1; i>=0; i--)); do run "${ARMS[$i]}" b; done
python3 $W/scripts/tsum.py $W/runs/$TAG
