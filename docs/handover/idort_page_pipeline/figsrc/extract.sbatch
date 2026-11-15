#!/bin/bash
#SBATCH -p general
#SBATCH -c 72
#SBATCH -t 01:30:00
#SBATCH -J idortfig
#SBATCH -o /raven/ptmp/jinma/idort_page/logs/extract_%j.out
export OMP_NUM_THREADS=1
cd /raven/ptmp/jinma/idort_page
/u/jinma/miniconda3/envs/magritte/bin/python /raven/ptmp/jinma/idort_page/figsrc/extract.py "$@"
