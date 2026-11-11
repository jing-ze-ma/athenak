#!/bin/bash
# run the athena binary ($ATH_BIN) as an MPI singleton detached from the slurm step's PMIx
# (otherwise one FATAL/MPI_Abort makes pmix kill the whole srun step and every other case)
for v in $(env | grep -o -E '^(PMIX|PMI|SLURM|OMPI|PRTE)[A-Za-z_0-9]*'); do unset $v; done
exec "$ATH_BIN" "$@"
