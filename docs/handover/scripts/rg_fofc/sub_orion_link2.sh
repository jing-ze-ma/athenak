#!/bin/bash -l
#SBATCH -J RG_fofc2
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
# Link 2 of the chain: restart from the LAST restart file job 1 wrote (nghost=3 is baked
# into it, so no override is needed and none would be honoured).
RST=$(ls -1 $D/rst/rg.*.rst 2>/dev/null | sort | tail -1)
if [ -z "$RST" ]; then
  echo "### RG_fofc chain: no restart file in $D/rst -- link 1 never got that far. Abort." >> out.txt
  exit 1
fi
echo "### RG_fofc chain link 2 restarting from $RST" >> out.txt
srun -n 32 --cpus-per-task=7 $D/bin/athena \
  -r $RST -i $D/rg.athinput -d . -t 03:55:00 >> out.txt 2>&1
