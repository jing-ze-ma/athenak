# NOTE caltech -> viper: Arepo2-IDORT GPU tests on Caltech H200, results (TASK ce65bc6e)

**Verdict.** Steps 1-3 pass. Step 4 stopped at a reproducible `nltree_mns_3d_cuda` CHECKS failure, so the nltree g25_red comparison was not run.
- **Not in git:** no Arepo source, diffs or patches were committed or pushed anywhere.
- **Full write-up (Caltech):** `/resnick/groups/carnegie_poc/jingze/arepo2_gpu/caltech/RESULTS_CALTECH.md`. The scripts and logs are next to it.
- **Where it ran:** the held allocation 5274909 (hpc-sm-02-04, 4 x H200, 32 cores), after run12 had left the node.

## 1. Builds
- **What was built:**
  - HYPRE 2.32.0 (CPU only);
  - the 10 Arepo configs from `gpu1/cfg` (idort-gpu-s1 6cef4d1d6);
  - 4 nltree builds (nltree-g25 1108335a71, in a separate clone).
- **Adaptations (only these):**
  - `Makefile.systype`: SYSTYPE stays Raven-OpenMPI, with three overrides:
    - CUDA arch sm_80 -> **sm_90**, in both the compile and the nvcc link line;
    - `-lbfd` dropped: there is no libbfd.so here, and it is only used by BACKTRACE, which no config sets;
    - `-march=native` -> `-march=znver4`: the build nodes are Intel, the H200 node is an EPYC 9554.
  - `-ffp-contract=off` and nvcc `-fmad=false` are untouched (in-tree).
  - **Modules:** gcc 13.2, cuda 12.9, spack gsl/fftw/hdf5 (serial), hpcx 2.17.1 Open MPI 4.1.7. A spack openmpi that the library modules pull in is unloaded first.
  - **Python:** a venv with numpy 1.26, h5py and scipy.
  - **Submodules:** the stale `lib/*/.git` gitlink files from submodules.tgz were removed. They point to a missing `.git/modules`.
  - **A7_rv:** the table paths point to the bundle tables.
  - **Gates:** they run in one srun step with SLURM_* cleared, so test.py's mpiexec sees the 32 cores as a plain host.
- **No defaults, preconditioner or solver settings were changed.**

## 2. Gates (IdortGpu 0 vs 1, same binary): all 10 ALL SAME
- **Cases:**
  - A1, A4, A4_m4, A6, A2;
  - I2a (1 island);
  - B3 (102/102 snapshots);
  - A7_rv;
  - I4_fhyb_levels_m4 (30/30);
  - I1_restart_bitwise.
- **What is identical:** snapshots, solver rows and budgets are bitwise equal, and check.py passes in both arms.
- **VERIFY lines are 0.** The verify call only exists on the IdortGpu 3 path; the GPU arms did run device solves.
- **I4_m4 matvec totals:** 4539 / 11368 / 6241 here, vs 4549 / 11349 / 6212 in Raven's older S0/S1a. This is cross-machine round-off, not a gate result.

## 3. g25_red timing
**Setup.** The arms ran interleaved in the order 0 1 1 0, with IdortGpuMinCells 1e5 and 4 islands x 8 = 32 ranks (Raven: 72).
- The solver rows are identical in all four arms.
- CPUsUsedPerGPU is 1, as in Raven's time_g1.sh.

| per full sync | Caltech H200 (32 ranks) | Raven A100 (72 ranks) |
|---|---|---|
| CPU arm | 8.69 s | 7.19 s |
| GPU arm | 4.17 s | 3.90 s |
| speedup | **2.08x** | 1.84x |

- **GPU arm, 2 full syncs:**
  - solve 6.33 s;
  - sweeps 2.44 s (1.22 s per full sync; Raven 1.39);
  - halo 0.51 s;
  - gpu io 4.16 s. This includes the init solve; Raven was ~2.3 s. The io profile was not run here.
- **44 partial steps:** about 8.2 s in both arms, on the CPU path (Raven 6.6 s).
- **The CPU arm is slower** because of the rank count (2.25x more cells per rank), not the hardware. A like-for-like CPU number needs a full 128-core allocation.

## 4. nltree-g25
**Passes:** `nltree_mns_3d` (CPU).

**Fails:** `nltree_mns_3d_cuda`, identically in 2 runs, at t = 1.334 of 2.0:
```
nltree_crosscheck_walk(), forcetree_nodelocal_checks.cc:176: NLTREE check: tree-walk acceleration of sink 0 deviates from the
brute-force reference beyond the opening-criterion bound at timebin 18: tree (-6.66418e-06, -2.86363e-05, 7.1693e-06) vs direct
(-1.10714e-05, -2.75734e-05, 2.96342e-07), |diff| 8.2335e-06 > 1.84603e-06 (gross 3.69207e-05)
```
- **Before the failure:** every GPU-vs-CPU parity check passed. The failure is the node-local walk against the brute-force bound on a tiny tree (14 sources, depth 2).
- **Logs (Caltech):** `arepo2_gpu/caltech/nltree/logs/test_nltree_mns_3d_cuda_n{2,3}.log`.
- **g25_red:** the CPU binary and the GPU binary (VORONOI_CELL_GPU + FORCETREE_NODELOCAL_GPU) are built, and `caltech/nl_g25.sh` is ready. That GPU build has no CHECKS. It was not run, because of the stop rule.

## Questions for viper / the user
1. Is the `nltree_mns_3d_cuda` crosscheck failure known on Raven, or is it new on H200?
2. Should the nltree g25_red comparison run anyway, without CHECKS?
3. Is a like-for-like CPU timing needed (a full-node allocation, which queues for days here)?

**Allocation.** 5274909 is free until about 10-25 06:12 PDT.
