#!/bin/bash -l
#SBATCH -J RG_fofc
#SBATCH -o /orion/ptmp/jinma/Athenak/red_giant/RG_fofc/log.out.%j
#SBATCH -e /orion/ptmp/jinma/Athenak/red_giant/RG_fofc/log.err.%j
#SBATCH --partition=p.exclusive
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=7
#SBATCH --time=04:00:00
module purge; module load gcc/13 openmpi/4.1
export OMP_NUM_THREADS=7 OMP_PROC_BIND=spread OMP_PLACES=cores
export UCX_LOG_LEVEL=error
D=/orion/ptmp/jinma/Athenak/red_giant/RG_fofc
cd $D
# FROM SCRATCH: FOFC+PLM needs nghost=3 and every red-giant restart bakes nghost=2.
srun -n 32 --cpus-per-task=7 $D/bin/athena \
  -i $D/rg.athinput -d . -t 03:55:00 >> out.txt 2>&1
