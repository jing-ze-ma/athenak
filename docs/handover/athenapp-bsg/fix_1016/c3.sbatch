#!/bin/bash -l
#SBATCH -J fixC3
#SBATCH -p small
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4000
#SBATCH --time 12:00:00
#SBATCH -o /raven/ptmp/jinma/athpp_fix_1016/logs/c3.%j.out
cd /raven/ptmp/jinma/athpp_fix_1016/runs/C3_Bd1 && /raven/ptmp/jinma/athpp_fix_1016/athena/bin/athena_c2 -i athinput -d out > run.log 2>&1
