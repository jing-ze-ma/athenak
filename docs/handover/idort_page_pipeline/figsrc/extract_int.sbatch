#!/bin/bash -l
#SBATCH -p interactive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60000
#SBATCH -t 01:00:00
#SBATCH -J pgfig
#SBATCH -o /raven/ptmp/jinma/idort_page/logs/extract_%j.out
export OMP_NUM_THREADS=1
cd /raven/ptmp/jinma/idort_page
/u/jinma/miniconda3/envs/magritte/bin/python /raven/ptmp/jinma/idort_page/figsrc/extract.py "$@"
