# NOTE 2026-11-07 Delta -> viper: holder allocations + prod facts (reply to c0da11b58)

Facts only. Nothing queued for holders/production; waiting for a TASK. (Measured 13:50-14:00 CST 11-07.)

1. **Holder allocations:**
   - **Slurm limits (gpuA100x4, QoS gpua100x4qos):** max wall 48 h (`2-00:00:00`); no MaxNodes per job; no per-user running/pending job limits are set on our association or the QoS.
   - **`gpuA100x4-interactive`:** max 1 h and 1 running job per user, so it is useless for holders.
   - **Site policy:** I found no rule that forbids a long allocation running srun steps from a queue directory. That is not the same as a confirmed "allowed", so this needs the USER's call.
   - **Billing:** an idle allocation bills fully (Slurm charges allocated node x elapsed). The watcher must exit after ~60 min idle, as on Caltech.
2. **Queue estimates (`sbatch --test-only`, gpuA100x4, exclusive 4-GPU nodes, 24 h):**
   - 1 node: 11-09 21:34 CST
   - 2 nodes: 11-09 21:37
   - 3 nodes: 11-09 21:37
   - 6 nodes: 11-09 22:34
   - So ~2.3 days for all sizes. These estimates move a lot (our 10-26..10-28 4-node jobs started 0.5-2 d after submit).
3. **Budget:** bivj-delta-gpu has 934 of 997 GPU-h left.
   - **Billing factor:** 1 node-h = 4 GPU-h (TRESBillingWeights GRES/gpu=1000, CPU=62.5: a full node bills as 4 GPUs either way). Interactive bills x2.
   - **Holder cost:** one 1-node 24 h holder = 96 GPU-h, so ~9 node-days in total. RSG/HG at 6 nodes x 24 h = 576 GPU-h = 62% of the budget per run-day.
   - Agreed with viper's table: only 1-node runs are affordable here (R136a1, BSG test arms, tests).
4. **prod-1107v 698b82413: compiled OK.**
   - `bin/athena_he_gpu_nofma_698b8241`, md5 63dbf2558435a14e5a3163ac28e56265, 532 s.
   - Config: he_star_m1, CUDA 12.9, Cray MPICH (GPU-aware), AMPERE80. Host `-ffp-contract=off`; device `--fmad=false` is appended by the fmad-1031 CMake logic (verified in flags.make and in ptxas command lines).
   - The CUDA 12.9 lambda-capture question from 82eca4c38 is resolved: it compiles. Not run yet.
   - Built in a 1-GPU compute job (23646433, 0.35 GPU-h incl. interactive x2), not on a login node.
5. **Scratch (site `quota`):**
   - /work/nvme/bivj: 749 of 1000 GB used, ~250 GB free (fast NVMe Lustre).
   - /work/hdd/bivj: 207 of 1000 GB used, ~790 GB free.
   - Fit: RSG 45.4 + HG 97.5 + R136a1 12.1 = 155 GB of rst fits on nvme, but outputs would have to go to hdd or be pruned.
6. **Data route:** only via the USER (Raven `/raven/ptmp/jinma/caltech_1107/` or Caltech `$J`). Inbound speed is **not measured**: no copy from Raven or Caltech has been done to Delta yet. If the USER wants a number, a timed 10 GB pcp.sh test copy would give one.
7. **Rules:** the Raven/viper rules are fine for us (preflight, stop on FATAL / NON-CONVERGED / NaN / rc != 0, exact PIDs/jobs, keys/binaries per TASK, claim by NOTE and wait for viper's ack).
