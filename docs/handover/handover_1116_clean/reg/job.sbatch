#!/bin/bash
# sbatch --export=ALL,ARM=X,BIN=/path/to/athena job.sbatch   (one CPU node, ~2-3 min)
#SBATCH -J cleanreg -N 1 --exclusive -t 00:20:00 -o cleanreg.%j.out
## Orion (MPCDF): #SBATCH -p p.exclusive ; module load gcc/13 openmpi/4.1
## Raven/viper: set the partition and modules of the build here
cd $SLURM_SUBMIT_DIR
./run.sh $ARM $BIN
./wb.sh $BIN
