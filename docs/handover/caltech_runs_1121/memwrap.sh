#!/bin/bash
O=$1; shift
if [ "${SLURM_LOCALID:-0}" = 0 ]; then h=$(hostname -s); nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits -lms 2000 | sed -u "s/^/$h /" >> $O/mem.log 2>/dev/null & fi
exec "$@"
