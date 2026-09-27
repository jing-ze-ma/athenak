---
name: kokkos-pin-4400-not-4602
description: DONE 2026-09-26 02:40 PDT -- kokkos pin bumped 4.4.00 (08ceff92) -> 4.6.02 (d8e9af03), pushed 2e028a38; all gates bitwise/round-off as expected; k46 ~4 % faster median cycle on H200
metadata:
  type: project
---

Found by viper's HIP gate (09-26, handover commit 4979b3a1).

- The pin was taken from `main`'s gitlink (08ceff92). Its commit message says "Updated Kokkos to 4.6.02", but the submodule's `CMakeLists.txt` says 4.4.0, and `git describe` gives 4.4.00.
- **Every Caltech binary up to 09-25 used Kokkos 4.4.00.** Viper uses 4.6.2 (bench/wt_rgbox/kokkos). Kokkos 4.4.00 has no `Kokkos_ARCH_AMD_GFX942_APU`, so viper cannot build against the pin.
- The tag 4.6.02 is already in the local clone (`/resnick/home/jingze/ATHENAK/athenak/kokkos`), so no download is needed.

**DONE 09-26 (2e028a38).** Gates:
- dhj CPU: bitwise.
- dhj GPU k44 vs k46, and 2 vs 1 GPUs: state bitwise; only the hst 3-mom cancelling column differs (~1e-3).
- He box CPU: identical. GPU vs CPU: same round-off as before.
- sph_atm: same as 4.4 (F2/F3 lateral fluxes are ~1e-12 noise, so "bitwise" means apart from those).
- sph_wedge: runs on CUDA with no nvcc fixes; round-off.

build_inc.sh now wipes the build dir when the kokkos commit changes. Binaries: builds/athena_k46{,bx,bi}_{cpu,gpu}.

**Original TODO (done)**
1. Once no agent is building (build_inc.sh symlinks the main checkout's kokkos), set the gitlink to 4.6.02: `git -C kokkos checkout 4.6.02`, then `git add kokkos` in the repo.
2. Rebuild fully: new build trees, or wipe `builds/inc_*`, because CMakeCache and Kokkos objects change.
3. Gates: dhj CPU old (4.4) vs new (4.6.02), expected bitwise; GPU smoke; GPU vs CPU values; 1 vs 2 GPUs.
4. Push, and tell viper it can then drop the symlink.

Lesson: verify a version from the source (version file or tag), not from a commit message.
