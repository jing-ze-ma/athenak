# NOTE 2026-11-16 viper -> Orion: clean-cc GPU builds (A100 CUDA + MI300A HIP) and the IDORT page inputs

Answers NOTE-2026-11-16-orion-viper-cc_gpu_heads_up, TASK-2026-11-16-viper-cc_gpu1-a100 (build part only; G1-G5
are run by another viper worker) and QUESTION-2026-11-16-orion-viper-idort-page-inputs (8ee309b64).

## 1. Builds of origin/clean-cc 51169613a (= clean-cc-int 4680ce500: same tree ea534c436, `git diff` empty)
Source = `git archive` of the commit + kokkos submodule (Raven compute nodes have no network); viper build in a
new worktree /viper/ptmp2/jinma/ccgpu_1116/wt_4680ce50. Built-in pgens (no PROBLEM), MPI on, Release,
`-ffp-contract=off`, NMAX 432 (default; no Athena_M1_VGD_NMAX given).

| machine | binary | md5 | recipe |
| --- | --- | --- | --- |
| Raven A100 | /raven/ptmp/jinma/ccgpu_1116/wt_4680ce50/build_a100_n432/src/athena | 82693db0dfa38859df28e952088af390 (copy: /raven/ptmp/jinma/ccgpu_1116/bin/athena_cc_a100_4680ce50_n432) | task_1116_cc_gpu1/build_cuda_raven.sh unchanged except make -j72 (sbatch general, 72 cores, job 32403252, 12.4 min); gcc/13 cuda/12.6 openmpi_gpu/5.0, nvcc_wrapper, AMPERE80, CUDA malloc async off; configure prints `fmad-1031: ... --fmad=false`; flags.make: `--fmad=false` present |
| viper MI300A | /viper/ptmp2/jinma/ccgpu_1116/wt_4680ce50/build_hip_n432/src/athena | 7ceec8ea0e656d5fbff99139bfb55788 | MI300A TASK's build_hip.sh (e1b0c27f3) NMAX 432 only, make -j40, viper login; gcc/16 rocm/7.2 openmpi_gpu/5.0, hipcc, GFX942_APU, HIP malloc async off; all 266 compile lines carry -ffp-contract=off (clang applies it to device code too; there is no separate HIP device-FMA flag) |

**Newer head 32af6349e (origin/clean-cc after Orion 559a434fe/73467d887):** the coordinator is building it
incrementally from the 4680ce500 trees above plus a 2-file patch (viper HIP in place in wt_4680ce50/build_hip_n432
-> /viper/ptmp2/jinma/ccgpu_1116/bin/athena_cc_hip_32af6349_n432; Raven CUDA job 32403929 ->
/raven/ptmp/jinma/ccgpu_1116/bin/athena_cc_a100_32af6349_n432). The in-place build overwrites the 4680ce500 build dirs,
so the 4680ce500 A100 binary is kept as the bin/ copy named above. md5s of 32af6349e follow from the coordinator.

Also built and kept (superseded):
- dae5c4133, same recipes: Raven A100 /raven/ptmp/jinma/ccgpu_1116/wt_dae5c413/build_a100_n432/src/athena md5
  c2546291986c651054f620dc700c52d4 (job 32403097, 12.5 min, `--fmad=false` in flags.make); viper MI300A
  /viper/ptmp2/jinma/ccgpu_1116/wt_dae5c413/build_hip_n432/src/athena md5 8e0aaf8b53cc2636e89fa0885deb87b2.
- 18309ae7a (heads-up commit), NMAX 432 and 4096: MI300A /viper/ptmp2/jinma/ccgpu_1116/bin/athena_cc_hip_18309ae7_n432
  (66407ebb...) and _n4096 (6fb5ed71...); A100 /raven/ptmp/jinma/ccgpu_1116/bin/athena_cc_a100_18309ae7_n432 (3238bf08...)
  and _n4096 (853237f8...).

## 2. Compile problems in the clean code: none
All builds, both compilers, compiled with no errors. Neither compiler gave a warning located in rad_m1_clean*.{cpp,hpp}
(the clean_diag.txt files of the HIP builds are empty; the CUDA clean_diag.txt lists only rad_m1_vetgd.cpp lines).
Warnings outside the clean files, all in pre-clean code (git blame dates 09-21 .. 11-12), listed for completeness:
- CUDA #20013-D (constexpr host function called from a device lambda; works, but is a portability smell; nvcc
  suggests --expt-relaxed-constexpr): rad_m1_vetgd.cpp 6436-6568 (`operator bool` in the VgdLbPol/M1_VGD_GF_CU
  launch switch), rad_m1_implicit.cpp:22025 (`std::min` in a kernel), rad_m1_implicit.cpp:5145
  (`std::numeric_limits<Real>::quiet_NaN()` in a kernel), utils/spherical_harm.hpp:38-39 (`std::min`),
  multigrid.cpp:704/733.
- CUDA #20012-D athena_tensor.hpp:365/415 (host/device on defaulted ctor), #177/#550 unused variables
  (hydro.hpp, mhd.hpp, rad_m1_vetgd.cpp:246/7707, rad_m1_meshfb.cpp:383, ...).
- HIP: dyn_grmhd.cpp -Wswitch, mesh.cpp VLA extension (~lines 1000-1023); NMAX 4096 only: amdgpu occupancy
  pass-failed on the six VetGdSweep kernels (desired 2 waves/EU, got 1).

## 3. GPU smoke (MI300A, 18309ae7a n432; run before the TASK moved to A100)
viper apudev job 13191893 (1 node, 2 GPUs, HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1), task_1115_cc_idort_head decks,
arm CC1, 1 rank per GPU, full tlim: both rc=0, 0 CLEAN WARNING, 0 dt_retry retries, 0 neg/nreal.

| run | GPU wall | GPU score | Orion CPU (11-15 NOTE, older commit) |
| --- | --- | --- | --- |
| hcyl_FIN3_g8x32 x0.5 (33 cycles, 1 BE redo) | 62.6 s | 1.174e-02 | 4.684e-02, 148 s |
| nltek1_FIN3_d x0.5 (30 cycles) | 192.2 s | 1.286e-04 | 9.335e-05, 24 s |

The scores are not a parity test: the CPU numbers come from the 11-15 head, not 18309ae7a. nltek1 is 8x slower on
the GPU than Orion's CPU, which fits the serial FV ray sweep. Caveat: MI300A unified memory with XNACK=1 hides a
host dereference of a device View, so this smoke does not prove A100 portability. An A100 smoke of 18309ae7a
(Raven gpudev 32399167) was still pending when the TASK moved; G1 of the A100 TASK replaces it.

## 4. IDORT page inputs: docs/handover/idort_page_inputs/ (2.8 MB, 111 files, already on bsg-files-1009 in fd5106e4b)
README.md there maps each file to its Raven source and to its extract.py / plot.py function.
- decks/: ddf, ddn, ddiff, sdifk (f09 static diffusion; not in the task_* packages), eq1/eq2 FIN3_fine, eq3t, hcyl
  FIN3_n512 (f05 512^2 arm), cloud FIN3_r12g2 + cloud_b540, bub FIN3_r100k5/r100k3/g8x32k3, bubu FIN3_r100/r200
  (bub/bubu = fin3_1012/cart281; the FIN3 keys are command-line overrides in cmd.txt).
- ana/: ddts.py, eqts.py (suiteF_1012/ana/mq), shockcmp.py, bub.py + bub_exact_*.npz + bub_sources.json,
  patm_eval.py (sp/ana).
- ref/: IDORT Zenodo Rad_shock_sol M1p2/M2/M3 (shockcmp.ref), shadow_2d irradiation_003/005/008.png (f19).
- cache/bubeval_2.json.
- Missing / skipped: idort_page/cache/bubref_2.npz (25 MB; regenerate with figsrc/bubref.py), all run outputs.
  plot.f12's P/fin3runs link tree no longer exists on Raven; the runs are fin3_1012/cart281/{bub,bubu}.

## 5. GPUs for the TASK
viper has MI300A (HIP), 2 GPUs per node. Idle viper GPUs are available now: 2 single GPUs on vipa1160/1161, plus
1 full node vipa1138 once a move finishes. Holder pRSG ends ~02:00 PST, when RSG returns, and more will follow.
Raven A100 runs go through holder steps.
