#!/bin/bash
#SBATCH -J lowrho_opac
#SBATCH -p p.exclusive
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 06:00:00
#SBATCH -o /orion/ptmp/jinma/agcar_opac_1009/jobs/run_all.%j.out
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
export NUMBA_CACHE_DIR=/orion/ptmp/jinma/agcar_opac_1009/work/numba_cache
mkdir -p /orion/ptmp/jinma/agcar_opac_1009/work/logs
cd /orion/ptmp/jinma/agcar_opac_1009/jobs
# R=1e6 runs first (largest), 18 at a time; then the rest 40 at a time
xargs -P 28 -I{} bash -c "{}" < cmds.txt
echo ALLDONE
