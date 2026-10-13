# NOTE 2026-10-13 from Caltech: fin3-cartlong -- NOT STARTED (no CPU slots), ready to go (re TASK-2026-10-13-caltech-fin3-cartlong.md)

- Binary built: fin3-1012 @ 281dd686, plain CPU (OpenMP + MPI, Release, default pgens) -> athena_cpu_omp_281dd686, md5 9ae9f09aeb81841a275e8d046a0eebe4.
- Files staged: /resnick/groups/carnegie_poc/jingze/cartlong_1013/caltech_cartlong (runs/ as pushed; run_cart.sh copy with /usr/bin/time -> a stand-in, not on Caltech compute nodes). Job: 1 node, 32 cores, 8 runs x 4 threads, 5 h.
- **Not submitted: no Caltech CPU capacity.** sbatch --test-only (10-12 18:00 PDT):
  - expansion (CPU) partition, normal or long QoS, any size (8-32 cores, 3-5 h): est. start **2026-10-21**;
  - debug QoS: max 30 min, 1 job -> too short (runs 0.5-2 h);
  - gpu partition, CPU-only on an L40S node (1 GPU idle, GPU-billed): est. start 2026-10-13 23:00 PDT.
- Proposal: leave these to Viper/DeltaAI (DeltaAI TASK 0f398db0 covers 4). Coordinator: reply here if you want the L40S-node run (or the 10-21 CPU queue) anyway; we submit within minutes.

# BSG FIN3 status (see NOTE-2026-10-13-caltech-bsg-fin3.md)
- On 281dd686 (md5 5cb31cf6); smoke + restart/cadence checks 4418238 pending, H200 est. 2026-10-13 11:20 PDT; links 4404881-83 behind it.
