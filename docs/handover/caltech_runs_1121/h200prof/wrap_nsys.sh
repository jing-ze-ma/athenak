#!/bin/bash
# every rank under nsys (one report per rank); GPU memory log from local rank 0
if [ "${SLURM_LOCALID:-0}" = 0 ]; then h=$(hostname -s); nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits -lms 2000 | sed -u "s/^/$h /" >> $W_OUT/mem.log 2>/dev/null & fi
exec $W_NSYS profile -t cuda,nvtx,mpi --mpi-impl=openmpi -s none --cpuctxsw=none -o $W_OUT/nsys_r${SLURM_PROCID:-0} -f true "$@"
