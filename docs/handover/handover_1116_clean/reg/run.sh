#!/bin/bash
# run.sh ARM BIN: the 29 runs of cmds.txt for the athena binary BIN in arm ARM (dirs
# named after ARM), NPAR at a time (default 30: one 72+ core node in ~2 min)
ARM=$1; BIN=$(readlink -f $2)
export REG=$(cd "$(dirname "$0")" && pwd) ARM BIN
export ATHENA_VIS=${ATHENA_VIS:-$REG/cc/ana}
export OMP_PROC_BIND=false OMPI_MCA_hwloc_base_binding_policy=none
cd $REG
mkdir -p rhse/$ARM && cp rhse/rhse_sp.athinput rhse/$ARM/
for c in uni closed thin1 feed1; do mkdir -p cs/${c}_$ARM && cp cs/$c/cs.athinput cs/${c}_$ARM/; done
mkdir -p cs/thin1mb_$ARM && cp cs/thin1/cs.athinput cs/thin1mb_$ARM/
for r in rw/N*_amp*; do [ -f $r/rw.athinput ] && mkdir -p ${r}_$ARM && cp $r/rw.athinput ${r}_$ARM/; done
SECONDS=0
xargs -d '\n' -P ${NPAR:-30} -I{} bash -c 'S=$SECONDS; {}; echo "END $((SECONDS-S))s :: ${1:0:120}"' _ {} < cmds.txt
echo "ALLDONE wall=${SECONDS}s"
