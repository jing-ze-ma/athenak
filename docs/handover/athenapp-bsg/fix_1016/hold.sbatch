#!/bin/bash -l
#SBATCH -J fixhold
#SBATCH -p general
#SBATCH -N 1
#SBATCH -c 72
#SBATCH --mem=120000
#SBATCH --time 24:00:00
#SBATCH -o /raven/ptmp/jinma/athpp_fix_1016/logs/hold.%j.out
sleep 86000
