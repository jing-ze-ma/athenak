#!/bin/bash
# rank 0 under ncu on ONE launch of the kernel whose mangled name is NCU_RE (exact), after NCU_SKIP matching launches; others plain
if [ "${SLURM_PROCID:-0}" = 0 ]; then
  exec $W_NCU --kernel-name-base mangled -k "$NCU_RE" --launch-skip $NCU_SKIP --launch-count 1 --set full \
       --import-source no -o $NCU_OUT -f "$@"
fi
exec "$@"
