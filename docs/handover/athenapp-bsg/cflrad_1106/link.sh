#!/bin/bash -l
#SBATCH -J repro3d_seeded
#SBATCH -A mpa_cpu
#SBATCH -p general
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=72
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH -o /raven/ptmp/jinma/athpp_repro3d/prod_seeded/link.%j.out
module purge; module load gcc/13 openmpi/5.0 hdf5-mpi/1.10.10
export LD_LIBRARY_PATH=$HDF5_HOME/lib:$LD_LIBRARY_PATH
cd /raven/ptmp/jinma/athpp_repro3d/prod_seeded || exit 1
[ -f STOP ] && { echo STOP file, exiting; exit 0; }
B=/raven/ptmp/jinma/athpp_repro3d/athena_exact
if [ -f output/bsg.final.rst ]; then
  cp -p output/bsg.final.rst output/bsg.last_link.rst
  srun -n 2304 --cpu-bind=cores $B -r output/bsg.final.rst -t 23:30:00 >> run.log 2>&1
else
  srun -n 2304 --cpu-bind=cores $B -i athinput -t 23:30:00 >> run.log 2>&1
fi
echo link done rc=$? $(date)
