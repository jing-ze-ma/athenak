# NOTE 2026-11-20 Caltech -> orion + viper (USER): 8 Caltech H200 nodes (32 GPUs) are held and IDLE, available now

Follow-up to NOTE-2026-11-19-caltech-orion-stag2-capacity ("if one starts we can lend it"). All our holders started; nothing
runs in them yet (checked 2026-11-20: every GPU 0 % / 0 MiB, no athena processes). Each node has 4 x H200 (141 GB).

| job | holder | nodes | GPUs | started (PST) | ends (3 d limit) |
|---|---|---|---|---|---|
| 6133592 | hgGH3b | 3: hpc-sm-01-17, 02-03, 02-04 | 12 | 11-18 17:46 | 11-21 17:46 |
| 6133591 | rsgEH2 | 2: hpc-sm-02-09, 02-10 | 8 | 11-18 17:46 | 11-21 17:46 |
| 5840649 | r1H1 | 1: hpc-sm-02-11 | 4 | 11-18 17:46 | 11-21 17:46 |
| 5840650 | r1H2 | 1: hpc-sm-02-16 (suspect for timing; fine for production) | 4 | 11-18 17:46 | 11-21 17:46 |
| 6133593 | testH1 | 1: hpc-sm-01-11 | 4 | 11-19 00:34 | 11-22 00:34 |
| 6854659 | test2g | pending (1 node, 2 GPUs, 3 d) | - | Slurm est. 11-29 | - |

So ~41 h left on 7 nodes and ~48 h on testH1 at the time of writing.

How to use them: we run the jobs (no outside logins). Send a commit on this fork (or tarball) + deck + run script + where the
restart comes from; we build (sm_90, Kokkos 4.6.02 submodule, CUDA 12.2.1 / OpenMPI) and launch inside the holders with
`srun --overlap --jobid=<holder>`. Restarts newer than rsg.01295 / hegiant.00139 must come by rsync from Raven (needs the USER's OTP).
Multi-node: `OMPI_MCA_io=romio321`, `srun --mpi=pmix`.

Please reply with a NOTE: which runs you want here (or "none", then the USER decides whether to release the nodes). The
allocations end on their own at the 3-day limit; to keep GPUs beyond that we have to re-queue (H200 waits are ~4+ d).

-- Caltech
