# NOTE 2026-11-16 Orion -> viper: gpufix regression (TASK d3ebb04ac)

## 1. clean-cc-a100fix-viper d6b6ed4e2: BITWISE to 62d71accc on every deck
- Full reg/ harness: all 29 cmds.txt runs, the 4 RT restarts (hcyl, sp3mb, feed, thin1; cycle-16 bins identical,
  "hydro_u same m1 same") and the 4 wb levels (8.39985e-08 6.83945e-08 5.56674e-08 1.85442e-07, = ref).
  cmp.sh differ 0 everywhere, in both directions.
- The CLEAN lines of all 37 logs (4014 lines) are byte-identical, and so are the wb logs. ana.sh = ref_62d71accc.txt.
- 0 CLEAN WARNING, all rc = 0.
- Binaries: gcc 13 / openmpi 4.1, Release, -ffp-contract=off, PROBLEM=built_in_pgens, VGD_NMAX 4096, kokkos 4.6.02.
  d6b6ed4e2 md5 e0899b8cb1789f4a5275c18492ee70ca, 62d71accc md5 a89bd2b77ccdf110ac528341c63ebd31.
  Job 230479 (2 nodes, both arms at the same time, 2 min). Dir /orion/ptmp/jinma/clean_1115/gpufix_reg.
- Harness notes: running the MPI binary without mpirun inside an srun step needs the SLURM_/PMI_/PMIX_/OMPI_ variables
  unset. The README's build line is missing -DPROBLEM=built_in_pgens. wb.sh always writes wb/runs, so copy the dir per arm.
=> CPU green for the push.

## 2. g3_rhse_sp3 uninitialised-read hunt: CPU is clean; valgrind not run
- **valgrind is not installed on Orion.** Building it from source was not allowed in this session. I used these
  substitutes on d6b6ed4e2, nlim 5 unless noted, comparing CLEAN it/sum/pic lines and both .hst by md5:
  - environment size: plain, env -i, +3000-char and +9000-char variables, on -O1 -g, Release and MPI+OMP builds
    (21 Release runs plus the -O1 set) -> all bitwise identical;
  - stack: -ftrivial-auto-var-init=pattern vs zero -> bitwise identical, and equal to the plain build;
  - heap: MALLOC_MMAP_MAX_=0 + MALLOC_PERTURB_ 165/90/1, serial, MPI 1 rank, and OMP 1 and 8 threads -> identical;
  - full run to tlim 2.5, perturbed (3000-char env + mmap off + perturb 165) vs plain -> all dumps bitwise;
  - -O2 -Wuninitialized -Wmaybe-uninitialized on rad_m1_clean*, rad_m1.cpp, vetgd, vetgd_cs, vetcol, vetlat, time2,
    tasks, update -> 0 warnings;
  - clang MSan: unusable (uninstrumented libstdc++ gives false reports).
  So no host uninitialised read in the CPU-compiled code affects this deck.
- **CPU head pic at cycle 1 = 3** (steady init 14 passes):
  `0/0:20)6 16)5 16)6 13)5 13)5 10)5 9)4 7)4 6)4 4)3 4)3 2)3 1)2 1)2 | 0/1:1)3 | 1/2:1)3 1/3:0)2 ...`
- **IMPORTANT: the task_1116_cc_gpu1 reference for g3_rhse_sp3 is NOT head.** Its ref/dumps/g3_rhse_sp3/rh.user.hst
  (md5 de8b1c18) is byte-identical to our 4680ce500 CPU run (commit 4e4345d67). 4680ce500 is 14 commits before
  62d71accc and lacks the Eq. 56-57 sources, the WB fix, torder and prod3. That ref has pic 6 at cycle 1 and 16 init
  passes. If your "pic 6 vs 3" or "1.2e-3" is measured against ref.json, the "6" mode may be the stale reference. Please
  check which A100 mode equals CPU head.
  CPU head reference to tlim 2.5: /orion/ptmp/jinma/clean_1115/gpufix_vg/runs/full_rel/ (bin/, rh.user.hst md5
  f94c1767, out.log). The env-size flip itself is still a real GPU-only nondeterminism.
- Where to look on the GPU side, since the CPU is clean (hypotheses):
  1. CUDA-only host paths the CPU never compiles: rad_m1_vetgd.cpp VgdClosureToDevice closures (~6372-6380, 6539-6560),
     VGD_CUDA_LB (~1280), the cub scan in VgdScanOne (~580); rad_m1_vetgd_cs.cpp shellf_/mainf_ (~2006-2115).
  2. Device uninitialised memory: `compute-sanitizer --tool initcheck` on A100, at both env sizes.
  3. If the A100 node has HMM (`nvidia-smi -q | grep -i Addressing`), a kernel reading host stack through a captured
     pointer/reference would depend on env size while the SASS stays identical. This is invisible on CPU. No obvious
     capture found in rad_m1_clean*.

-- Orion (CPU help only; viper owns origin/clean-cc)
