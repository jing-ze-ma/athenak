# TASK for the DeltaAI session: bring up AthenaK on NCSA DeltaAI (2026-09-28)

You are a Claude Code session on **DeltaAI** (NCSA). The user develops AthenaK on two other machines:
- **viper** (MPCDF, AMD MI300A APUs, HIP);
- **Caltech Resnick HPC** (NVIDIA H100/H200, CUDA, x86 hosts).

DeltaAI is a third GPU platform: NVIDIA **GH200 Grace Hopper** nodes, with an ARM (aarch64) Grace CPU and an
H100-class GPU per superchip. The user wants to know whether AthenaK builds, runs correctly and performs well
there. This file is your whole brief. Read it to the end before doing anything.

## 0. Getting the code (no existing checkout needed)

You do not need the full history. Clone ONE branch, shallow, with the Kokkos submodule:

```bash
git clone --depth 1 --branch rt-integration --single-branch --recurse-submodules --shallow-submodules \
    https://github.com/jing-ze-ma/athenak.git athenak
```

- In-source builds are refused by CMake.
- Kokkos is the in-tree submodule (4.6.02). Do not use a site Kokkos.
- Data (correlated-k tables, He-box files) lives on a separate data-only branch:
  `curl -LO https://github.com/jing-ze-ma/athenak/raw/data-2026-09-26/athenak_data_2026-09-26.tar.gz`
  (md5 057e168a9a1559e9eb419bf362dae0b1). Unpack it outside the repo. The inputs use the placeholders
  `ATHENAK_CK_DATA` (the `exo_fms_ck/` dir), `ATHENAK_CK_DATA10` and `HE_BOX_DATA`. Substitute the real paths
  with `sed` into a copy of the input; never edit the tracked inputs.
- Read `CLAUDE.md` at the repo root (build, test, style, architecture).
- Do the one-time credential cleanup in `docs/handover/NOTE-2026-09-28-token-cleanup.md` first.
- Pushing needs GitHub credentials. Until the user sets up an SSH key here, do not push. Write results to files
  and tell the user; the user relays them.

## 1. Find out the machine first (verify, do not assume)

Everything below about DeltaAI is from general knowledge and must be checked on the machine: `module avail`,
`sinfo -s`, `sacctmgr show assoc user=$USER format=account,partition,qos`, `nvidia-smi` on a compute node,
`uname -m`, and the NCSA DeltaAI user docs.
- **Nodes:** expected 4 x GH200 per node (aarch64 Grace + sm_90 GPU, coherent CPU-GPU memory).
- **Batch:** Slurm. Find the GPU partition names (e.g. `ghx4`, an interactive/debug one), the account string
  and the per-job limits.
- **Toolchain:** the CUDA toolkit module, the MPI that is CUDA-aware (probably Cray MPICH with
  `MPICH_GPU_SUPPORT_ENABLED=1`, or an OpenMPI/HPC-X build), and a gcc that nvcc accepts.
- **Launcher:** plain `srun` vs `srun --mpi=pmix`. On Caltech the default pmi2 started singleton ranks.
  Check that 2 ranks really see each other: the log prints the rank count.

Record what you find in `DELTAAI_FACTS.md` next to your work dir. It is the most useful single output of the
first session.

## 2. Build

The GPU build recipe to adapt is `docs/handover/caltech-2026-09-26/scripts/build_caltech.sh`. It builds from
a `git archive` snapshot, so a dirty tree never leaks into a binary. Use the same Kokkos flags plus the Grace
host arch:

```bash
cmake -B build -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D Kokkos_ARCH_ARMV9_GRACE=On \
  -D CMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper \
  -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D PROBLEM=deep_hot_jupiter_rt
```

- Also build `-D PROBLEM=box_convection` (He box, M1 radiation) and one with no PROBLEM (built-in test pgens).
- If `Kokkos_ARCH_ARMV9_GRACE` fails with the host compiler, drop it and note it. It only affects host code.
- Build on a compute node if the login nodes limit memory or cores. Keep build directories outside the repo.
- Record the modules, the compiler versions and the md5 of each binary.

## 3. What to run, in this order (report after each step)

1. **Smoke:** `./athena -c` (build configuration), then one tiny run per binary (a built-in 3-D linear wave;
   see `tst/test_suite/`), on 1 GPU.
2. **Test suite:**
   ```bash
   cd tst
   python run_test_suite.py --gpu "-DKokkos_ARCH_HOPPER90=On -DKokkos_ARCH_ARMV9_GRACE=On -DCMAKE_CXX_COMPILER=$PWD/../kokkos/bin/nvcc_wrapper"
   ```
   Run it inside a GPU job. Report pass/fail per test. A failure here matters more than any timing.
3. **Rank-count consistency:** WASP-121b 1x input (`docs/handover/caltech-2026-09-26/inputs/wasp121_1x/`,
   `sparc_w121_1x.athinput` + `ic_w121_1x.txt`, ck paths substituted). Run about 50 cycles on 1 GPU and on
   2 GPUs from the same fresh start. Compare **rst/bin payloads**, not hst bytes: the hst 3-moment column is
   a cancelling sum that differs at about 1e-3 between rank counts. Expect bitwise equality.
4. **Timing**, same binary, all arms interleaved in ONE job; report the median of 8-cycle windows:
   - The production-grid variant (nx1 = 256, the `G_PROD` line in `grid_w121_1x.env`), 1 GPU and 2 GPUs,
     production keys as in the input.
   - Caltech H200 reference for this benchmark: **21.56 ms/cycle on 2 GPUs, 36.37 ms/cycle on 1 GPU**
     (rt-integration 8be0ad67).
   - Also try 4 GPUs (one node) if the mesh decomposes (24 MeshBlocks on the cubed sphere; <= 24 ranks).
5. **He box M1 (optional, if 1-4 are clean):** `box_convection` with the He-box files from the data tarball.
   Take the input from `/viper` only via the user, since it is not in the repo. Otherwise skip it.

## 4. Known state of the code (important)

- **Do NOT start any WASP-121b / dhj production run.** A bug is being fixed on viper right now: with
  `rot_potential = true` and `etotgrav = true`, the horizontal centrifugal work is counted twice in the energy
  equation (`src/pgen/deep_hot_jupiter_rt.cpp` SourceFunc), which heats the deep atmosphere. Short tests,
  gates and timings are unaffected in any way that matters for this bring-up.
- Other fixes are in flight on viper branches (He-box wall potential, FOFC with etotgrav, energy guards).
  Test `rt-integration` as it is. Do not try to fix physics here.
- GPU-bitwise notes: Caltech made one CUDA rounding explicit (`CkLayW`, `__fma_rn` under `__CUDA_ARCH__`) so
  that nvcc keeps H100/H200 results bitwise across the ck-lin2 refactor. GH200 is sm_90 too, so a bitwise
  mismatch between rank counts on DeltaAI would be a real finding. Report it with the first differing
  variable and cell.

## 5. Working rules (from the user; same as on the other machines)

- **Agents:** Opus 5.5 at medium effort for delegated work (no haiku/sonnet). One agent per branch/worktree.
  Agents do not merge or push.
- **No blocking waits:** do not run `sleep`/poll loops in the foreground. Submit, look once, return control.
  One background watcher per batch of jobs, polling `squeue` every 10 minutes.
- **Save tokens:** narrow tasks, no surveys, one deliverable per agent.
- **Timing comparisons:** same binary, interleaved arms, one job. Correctness first, speed second; never trade
  accuracy for speed.
- **Tight time limits** on batch jobs, so backfill can start them.
- **Never write into `run/`** in the repo. Never commit data or binaries.
- **Style:** CI enforces cpplint with 90 columns (see `CLAUDE.md`), if you ever touch code.
- **Credentials:** never put tokens in remotes, files or commits. Use an SSH key or `gh auth`.

## 6. Deliverable for the user

A short report, also written to `DELTAAI_REPORT.md` in your work dir:
- the machine facts;
- build recipe and status per problem;
- tst results;
- 1-vs-2-GPU bitwise result;
- ms/cycle for 1/2/4 GPUs next to the H200 numbers above;
- anything that needed a workaround.

If the user sets up push access, commit `DELTAAI_FACTS.md`, `DELTAAI_REPORT.md` and your build/run scripts
under `docs/handover/deltaai-2026-09-28/` on a NEW branch `deltaai-bringup` (never on rt-integration). The
viper session merges it after review.
