# NOTE 2026-11-17 Orion->viper (USER): race audit part 2 (rad_m1, other kernels, stag2 scratch table)

This follows ab09c9e2b and answers d001a877f. It was done by code reading.
- **Trees:** clean-cc-stag2-gpu2 (6c138e42b), gpu1 (9689f4216), clean-cc (dfe21a746), prod-1118 (b2d852bb6), rsg-ck-1008 (19fb0b749).
- **Independent agreement:** we found R1 (hydro and MHD scalar WAR), R2 and R4 on our own, in all trees (rsg-ck-1008 included).
- **Our MHD patch:** +2 lines, the same placement as yours.

## (a) Results not already in d001a877f
- **rad_m1: no findings in any tree.**
  - The staggered branch changes (clean.cpp, clean_solver, time2, dtretry_visit, rad_m1.cpp) contain no team code.
  - gpu1/gpu2 split clean_cell into CleanCellLineSolveRows and CleanCellLineFactorSolve: persistent teams, rows in a private
    global cln_lw slice, PCR + Thomas.
    - Barriers sit after the load, after each PCR round, after Thomas and after the output.
    - The barrier after the output is the one needed for the write-after-read hazard on the next line's load. It is present.
  - CleanFvLevels is your C1 (benign).
  - The shared implicit (PCR, M1PCRX, stot/stopt, Det/LayReduce), precond, vet, vetcol, vetgd, vetgd_cs and launch kernels are
    all clean.
- **R3:** our sub-auditor did NOT flag it. On re-check we agree: formally a barrier is missing, and in practice the preceding team
  parallel_reduce syncs on CUDA/HIP.
- **O2, the same class as R3 (theoretical):**
  - Kernel: ck_impl_fused in two_stream_rt_driver.cpp.
  - With t0rec_, t0_(m,k,j,i) is written at clean 278 / prod 305 (rsg column_ck.hpp 1508). The row pass reads it at i±1 (619 / 661 /
    1864).
  - The only sync between them is the sg/ss team parallel_reduce. An explicit barrier exists only with ck_impl_cvsec or
    ck_impl_glob.
  - Reached with ck_impl_fuse + ck_impl_reuse_jac = 2.
  - Optional fix: +1 `tm.team_barrier()` after `}, ss);` (302 / 344 / 1547).
- **O1 (confirmed, not built in production):**
  - Where: pgen/deep_hot_jupiter_rt_old.cpp :846, :1042, :1896 (prod/rsg).
  - What: a serial column body in par_for_outer with no inner loop, so every lane runs it and `u0 += ...` (:987, :1193, :2116) is
    applied T times without atomics.
  - Reach: only -DPROBLEM=deep_hot_jupiter_rt_old. We found no such build.
  - Fix: wrap the body in Kokkos::single.
- **Everything else in scope is clean:** two_stream partition / beam / cad, red_giant, update / WB, coordinates, bvals, flux_correct,
  prolongation, diffusion and refinement_criteria. Shearing box has only your E2.

## (b) stag2 team scratch vs N (nx1 = interior cells, 8-byte Real)
- **Staggered clean kernels: 0 B.** CleanCellLineSolveRows, CleanCellLineFactorSolve and CleanFvLevels use global memory.
- **Other kernels on the stag2 clean path:**
  - vet_gd moments (m1_vgd_mom, vgd_hr, vgd_hrcs, vgdcs_mom): 8 KB, fixed.
  - vet_gd up_sum: 2.3 KB.
  - vet_col m1_vcol_team: 8(5 nx1 + (nrs+lc) nray + 24 lc) B, with lc capped so the total stays at about 3072 reals. That is
    ~4 KB at nx1 = 64 and ~24-26 KB at nx1 = 548.
  - **All level 0** on MI300A (~64 KB) and on A100 (48 KB).
- **Shared implicit M1, if that path is used:**
  - PCR (ImplicitPCRSolve, M1PCRX, launch): 64 nx1 B, or 80 nx1 cyclic. That is 35 / 44 KB at nx1 = 548.
  - M1PCRMG: 64 nx1 B. m1_gf_fac: 60 n1 B, which is 33 KB.
  - These are level 0 on MI300A and below 48 KB on A100, but above the 32 KB probe threshold of TeamScratchLevel.
- **Only level-1 case:** vet_mb VetMomOut, 128 nx1 B. That is 70 KB at nx1 = 548, so level 1 on both GPUs (A100 for nx1 >~ 385,
  MI300A for nx1 >~ 510). It is not on the clean path.

## Ownership
- viper owns the clean-cc fixes; Orion picks them up at rebase.
- The prod-1118 / rsg-ck-1008 patches (hydro +5, MHD +2) are held on Orion for the USER to decide.
- No code branch was touched.
