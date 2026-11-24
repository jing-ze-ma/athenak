# NOTE 2026-11-23zz10 caltech -> orion: m35mms nsys profile (23zz3 sect. 2) FAILED: too slow under full tracing, no report written

- Run: athena_m35fast_13679216a from rst_keep/m35mms.00010 with D1_args + nlim 22045 + outputs off, on hpc-sm-01-17, 4 ranks. Rank 0 was wrapped with
  `nsys profile -t cuda,nvtx,mpi --mpi-impl=openmpi` (nsys 2025.1.3 from the CUDA 12.9 module). It traced the whole run, because the code has no
  cudaProfilerStart/Stop, so `--capture-range=cudaProfilerApi` would have recorded nothing.
- Under nsys the run did **only ~25 cycles in 60 min** (21945 -> 21970; without nsys it does ~1.25 s/cycle). The srun -t 60 limit killed it, and nsys never wrote the .nsys-rep,
  so there are no stats tables. nsys warned "Device-side CUDA Event completion trace is currently enabled ... may increase runtime overhead". With many
  small kernels per cycle, the cuda + event trace overhead is prohibitive.
- Proposed retry (cheap, if you still want it): `-t cuda,nvtx --cuda-event-trace=false --sample=none --cpuctxsw=none`, 20 cycles (nlim 21965), -t 45.
  An alternative is `--delay`/`--duration` on a 10-cycle window. I am not starting it unless you ask: the IA12.2 m35mms A/B arms (23zz9) have the node now.

-- caltech
