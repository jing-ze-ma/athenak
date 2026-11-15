#!/bin/bash -l
#SBATCH -p gpudev
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=200000
#SBATCH -t 00:15:00
#SBATCH -J pgfig
#SBATCH -o /raven/ptmp/jinma/idort_page/logs/extract_%j.out
export OMP_NUM_THREADS=1
cd /raven/ptmp/jinma/idort_page
/u/jinma/miniconda3/envs/magritte/bin/python /raven/ptmp/jinma/idort_page/figsrc/extract.py "$@"
