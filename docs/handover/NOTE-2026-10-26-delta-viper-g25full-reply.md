# NOTE 2026-10-26 Delta -> viper (USER): reply to 43f11c08: Delta CAN run g25_full on 4 gpuA100x4 nodes; ~20-30 min of login-node build first (x86 GMP + HYPRE); nothing started

**Short answer.**
- **Yes:** 4 gpuA100x4 nodes (16 x A100-40GB, 256 cores), the 3 arms in one <= 1 h job.
- **Not blocked.** The x86 GMP and HYPRE builds must be done first (login node, no GPU-h).
- **Bundle:** verifies here: md5 ee4103cf matches, and `git bundle verify` against DeltaAI's arepo2n is read only and passes. The prerequisites 1108335a7 and 6cef4d1d6 are present.
- **Nothing started**, as asked.
- **DeltaAI already has a GO** for the same measurement (af4d2971). Delta adds value only as a backup, or as an A100 point with the same node layout as Raven.

## 1. Build (x86 + sm_80)
- **Modules present:**
  - cuda/12.9;
  - cray-mpich 9.1.0;
  - cray-hdf5 / cray-hdf5-parallel 1.14.3.9;
  - gsl/2.8-gcc13.3.1;
  - cray-fftw 3.3.10.11.
- **Not present as modules:** GMP and HYPRE. DeltaAI built both locally, but for aarch64, so they cannot be reused.
- **Plan:** rebuild both from the tarballs already on /work:
  - gmp-6.3.0.tar.xz;
  - hypre-2.32.0 with DeltaAI's build_hypre.sh.
- **Where:** into `delta_1008/arepo2/`. Clone from DeltaAI's arepo2n plus the bundle, and only read DeltaAI's tree.
- **Time:** ~20-30 min on the login node for GMP, HYPRE and the 3 Arepo builds plus the B3 build.
- **Flags:** `-arch=compute_80 -code=sm_80`.
- **Known Delta pitfall:** the Cray CC wrapper adds -march=znver3, and FMA contraction changed AthenaK results here. This matters only for bitwise comparisons, not for timing. If you want identical flags across arms, say so.

## 2. Queue and cost (sbatch --test-only, 10-26 11:04 CDT; 4 nodes x 64 tasks, 1 h)
| partition | est. start | GPU-h charged (1 h) |
|---|---|---|
| gpuA100x4-interactive (max 1 h, 4 nodes) | 11:16 (~12 min) | 32 (2x rate) |
| gpuA100x4 (regular) | 12:22 (~1.3 h) | 16 |

- **Queues now:**
  - regular gpuA100x4: 1710 pending / 232 running;
  - interactive: 9 pending / 20 running.
- **Fairshare:** 0.29 (was 0.047 on 10-09).
- **Caveat:** the --test-only estimates have been optimistic for 4 nodes before. On 10-09, 4-node interactive jobs waited up to ~1.5 d.
- **Balance:** **96 GPU-h** (not 99). Either partition fits, with >= 5 GPU-h kept in reserve.
- **Delta's recommendation:** use the regular partition if the ~1 h wait is acceptable. It costs half the GPU-h.

## 3. Wall time and memory (estimated)
- **Per node, Delta vs Raven:**
  - GPUs: the same, 4 x A100-40GB;
  - cores: 64 (Milan) vs 72;
  - host memory: 251 GB usable.
- **Ranks:** at 64 ranks per node, the CPU arms have ~11% fewer cores than Raven 4 nodes.
- **Timing:**
  - per top step, expect about the Raven 4-node numbers: GPU ~21 s, CPU ~28-32 s;
  - per arm, with init and 5 top steps: ~4-6 min;
  - 3 arms: ~15-25 min, which fits 1 h.
- **Memory:**
  - host: Raven peaked at 135 GB per node at 72 ranks, which fits 251 GB;
  - GPU: ~29 GB per GPU of 40 GB.
  - MaxMemSize must be <= ~3500 MB at 64 ranks per node. Please give the value Raven's run_full.sh uses.
- **Not needed on Delta:** the MPICH_SHARED_MEM_COLL_OPT hang was DeltaAI-specific. Delta would use plain `srun -K -n N --cpu-bind=cores`, because `--gpus-per-task/--gpu-bind=closest` hangs Cray MPICH GPU IPC on Delta (NOTE fb16dab9).

## 4. Needed before a start
- The user's site choice.
- The run_full.sh param and the build flags (comb/build2.sh). DeltaAI's copy is fine; say where it is.
- No Arepo source in git.
