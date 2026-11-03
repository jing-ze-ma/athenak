#!/bin/bash -l
# run_li.sh <workdir> [li.py args...]: load the CPU MPI environment and run the suite
# (Caltech: launch ON a compute node of a held allocation, detached:
#   ssh <node> "setsid nohup bash -l <this> <work> --only ... > <work>.out 2>&1 < /dev/null &")
module purge >/dev/null 2>&1
module load gcc/13.2.0-gcc-13.2.0-w55nxkl openmpi/5.0.1-gcc-13.2.0-xpoh5uw >/dev/null 2>&1
H=$(dirname $(readlink -f $0)); W=$1; shift
export OMPI_MCA_rmaps_base_oversubscribe=1 PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe
exec python3 $H/li.py --bins ${LI_BINS:-$H/bins.caltech.cfg} --work $W "$@"
