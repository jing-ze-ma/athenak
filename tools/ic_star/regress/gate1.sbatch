#!/bin/bash -l
#SBATCH -J icsgate
#SBATCH -p gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=18
#SBATCH --mem=120000
#SBATCH --time=00:59:00
#SBATCH -o /raven/ptmp/jinma/icstd_1016/logs/gate.%j.out
# ic-standard-1015 seam task: gate (a), 1 A100, 1 rank, seed 0, 3 d (as bsg_newsetup run1g.sbatch). env NAME, DECK, KEYS
export SBATCH_EXPORT=NONE KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1
module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0
W=/raven/ptmp/jinma/icstd_1016; BN=/raven/ptmp/jinma/overlap_1015/bin/athena_he_a100_b16b66db
R=$W/runs/${NAME:?}; mkdir -p $R; cd $R
echo "$(md5sum $BN) -i ${DECK:?} $KEYS" > cmd.txt; cat cmd.txt; hostname; nvidia-smi -L
srun -n 1 --cpu-bind=cores --export=ALL $BN -i $DECK $KEYS -d $R < /dev/null > run.log 2>&1
echo "rc=$? $(grep -c FATAL run.log) FATAL; $(grep 'cycle=' run.log | tail -1 | cut -c1-80)"
