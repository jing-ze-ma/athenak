#!/bin/bash -l
#SBATCH -J col4x4
#SBATCH -A mpa_cpu
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH -o /raven/ptmp/jinma/athpp_repro3d/col4x4/link.%j.out
module purge; module load gcc/13 openmpi/5.0 hdf5-mpi/1.10.10
export LD_LIBRARY_PATH=$HDF5_HOME/lib:$LD_LIBRARY_PATH
cd /raven/ptmp/jinma/athpp_repro3d/col4x4 || exit 1
srun -n 32 --cpu-bind=cores ../athena_col4x4 -i athinput -t 05:45:00 >> run.log 2>&1
echo done rc=$? $(date)
