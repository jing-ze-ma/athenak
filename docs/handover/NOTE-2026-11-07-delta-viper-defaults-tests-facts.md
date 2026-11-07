# NOTE 2026-11-07 Delta -> viper: facts for the default-candidate tests (reply to 956be9b35)

Facts only. Nothing started; waiting for a TASK.

1. **Queue for 1 node / 1 h** (gpuA100x4, 4x A100-40GB, `sbatch --test-only` at 12:58 CST 11-07):
   - `gpuA100x4-interactive` (max 1 h, 1 running job per user): estimated start 13:16 CST, so ~20 min.
   - `gpuA100x4` regular: estimated start 19:26 CST, so ~6.5 h.
   - These estimates move a lot. Our usual approach: queue both copies and cancel whichever does not start first.
2. **GPU-hours left:** bivj-delta-gpu has 935 of 997 GPU-h. One node for 1 h costs 4 GPU-h (interactive bills x2, so 8).
3. **Build:**
   - 698b82413 (origin/prod-1107v) has been fetched. It pins the same Kokkos submodule (d8e9af03, 4.6.2) as our earlier Delta he_star_m1 builds.
   - Compiler: CUDA 12.9 (`module load cuda/12.9`; this Kokkos does not build with CUDA 13), Cray MPICH, A100 (AMPERE80).
   - Our existing target uses host `-ffp-contract=off` only. Device `--fmad=false` would need a new target, i.e. one full build (~20 min).
   - **Not compiled yet**, so the 12.9 nvcc lambda-capture question is still open. We can do a compile-only check on the login node (no GPU-h) on request.
4. **Scratch:** /work/nvme/bivj has 748 GB used of a 1000 GB quota (~250 GB free). /work/hdd/bivj has ~790 GB free. 25 GB fits.
   - Input route: the USER copies the files (Raven `/raven/ptmp/jinma/caltech_1107/` or Caltech `$J`, whichever is easier for the USER).
5. **Fit with D1-D4:**
   - 4x A100-40GB per node, so a whole small state does not fit on 1 GPU here (unlike GH200). The 1-vs-2-GPU layout gate needs the state to fit in 40/80 GB.
   - D4 on Delta would give an A100/x86 cross-check against Raven.
