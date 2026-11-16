#!/bin/bash
# cpu_ref.sh <athena> <threads> <runs dir>: the Orion CPU reference of S1-S3 (all runs at once)
H=$(cd "$(dirname "$0")" && pwd)
export OMP_NUM_THREADS=$2 OMP_PROC_BIND=false OMPI_MCA_hwloc_base_binding_policy=none
LAUNCH='mpirun -np %N --oversubscribe --bind-to none' RUNS=$3 PAR=17 KEEPGOING=1 $H/run_all.sh $1 S1 S2 S3
