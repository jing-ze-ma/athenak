#!/bin/bash
# nsyswrap.sh OUTROOT BIN ARGS...: rank 0 under nsys (cuda, nvtx, mpi; Kokkos labels need KOKKOS_TOOLS_LIBS=<libkp_nvtx_connector.so> if available), others plain + memwrap
O=$1; shift
if [ "${SLURM_PROCID:-0}" = 0 ]; then
  h=$(hostname -s); nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits -lms 2000 | sed -u "s/^/$h /" >> $O/mem.log 2>/dev/null &
  exec nsys profile -t cuda,nvtx,mpi --mpi-impl=openmpi -s none --cpuctxsw=none -o $PWD/nsys_r0 -f true "$@"
fi
exec $(dirname "$(readlink -f "$0")")/memwrap.sh $O "$@"
