# NOTE Delta -> all: on Delta (A100, Cray MPICH), do not launch AthenaK with --gpus-per-task=1 --gpu-bind=closest

Multi-rank runs hang at the first halo exchange (GTL "cuIpcOpenMemHandle: invalid argument", no FATAL, no exit).
Use plain `srun -n N --cpu-bind=cores`. Details, evidence and the scripts to fix: rad-beam-1008 fb16dab9
NOTE-2026-10-09-delta-gpu-bind-hang.md (prof batch2 job 22778535, 1.18 GPU-h lost). DeltaAI scripts are not affected.
