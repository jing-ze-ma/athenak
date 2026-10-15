# NOTE 2026-10-09 (viper): production line merged into rt-integration

## What merged
- Merge commit df5be1ce: `fork/mem-1009` 6c5d8fb2 into `rt-integration` (ba7212d6), non-ff, no conflicts.
  It brings rc-1009 (36f18f1d), accel-1009 (incl. cherry-picks a07cf032 / b93b6a9d), hrup-1009 (half-range VET),
  truerepro2-1009 (BSG pgen keys), hrup-bsg-1009 (98835d99 rank-ordered shell sums, 26994b9a halo wait),
  and vet_gd_twin_lowmem (6c5d8fb2).
- After the merge, `src/` equals 6c5d8fb2's `src/` exactly. So the merged he_star_m1 binary has the same code as the
  Caltech BSG binary. The viper ROCm 7.2 build of the merge is also byte-identical to the 6c5d8fb2 build
  (md5 516194d3).
- 28dd7b0e: style only. One `}}` in rad_m1_implicit.cpp (from 83acafb4, 10-02) is split onto two lines, so the code
  does not change.
- No defaults changed in this step.

## Gates (all on the combined code)
1. Style. The repo-wide check fails on the old code as well (cpplint master, many old files). For the 26 src files
   the merge touched, cpplint counts per file are equal before and after the merge, so the merge adds no new lint. The
   custom checks (tabs, `}}`, `#pragma`, mode 644) find only the `}}` above, which is now fixed.
2. CPU regression, login node, 16 cores, logs in /viper/ptmp2/jinma/merge_1009/logs/regress_merge.log:
   - Tests run: rad_m1, hydro, mhd and rad (the GR family is skipped).
   - `_cpu`: 62 passed, 3 skipped.
   - `_mpicpu`: 4 passed, 3 skipped, 1 failed. The failure is rad/test_rad_lwave2d_amr_mpicpu, a segfault (rc 139).
     It fails the same way on the old rt-integration ba7212d6 (regress1_old_lwave2d.log), so the merge did not cause it.
3. Keys-off bitwise against the old rt-integration binaries (ba7212d6). CPU binaries, files in
   /viper/ptmp2/jinma/merge_1009/cpugate:
   - WASP-121b cubed-sphere MHD input (w121_mhd_1003/cpugate/off.athinput, 12 ranks, 4 cycles): hst byte-identical,
     dump data bitwise identical in 6 of 6 files.
   - He box M1 input (bench_hebox_0929, 4 ranks, 3 cycles): hst byte-identical, dump data bitwise identical in 5 of 5
     files.
   - The bin headers differ only by the parameter dump of newly recorded defaulted keys (sp_x2_periodic_image,
     implicit_flux_faces, ...).
4. Production bitwise against 98835d99, on Raven gpudev with 4 A100 (job 31033897, build 31033896; merged binary
   /raven/ptmp/jinma/merge_1009/bin/athena_he_a100_df5be1ce_merge1009, md5 452f281e). Both runs passed:
   - AG Car A hr production input, 10 cycles: 26 of 26 files (bin + hst) byte-identical.
   - BSG true-repro hr production input, reduced 4-block mesh, 5 cycles: 30 of 30 files byte-identical.

## Binaries
- viper: /viper/ptmp2/jinma/builds/bin/athena_he_gpu72_df5be1ce_merge1009 (= athena_he_gpu72_6c5d8fb2_mem2),
  athena_dhj_cpu_df5be1ce_merge1009, athena_box_cpu_df5be1ce_merge1009.
- Raven: /raven/ptmp/jinma/merge_1009/bin/athena_he_a100_df5be1ce_merge1009.
- Production stays as it is: Caltech BSG runs 6c5d8fb2 (same src), and AG Car / He giant run 98835d99. A swap to the
  merge does not change results (gate 4).

## Not merged (side branches)
- fluxmean-1009 a7ed7b3a (local only, not on the fork): `flux_opacity = rosseland | fluxmean` force option, default
  rosseland. It is open work for the AG Car flux-mean table and is not superseded.
- fsanchor-1009 2edf2d23: implicit_fs_anchor (WIP), implicit_fs_lateral and vet_gd_fs_moments. It also carries the
  rad-beam commits below. Experimental and still open.
- rad-beam-1008 c7342c4d (fork tip is now 1f582486): code commits 29a47fae (vet_col_lat_fallback = local) and
  382389a6 (implicit_blend_ffs), not in rt-integration; the rest is Delta docs. vet_source_noesrc from this line is
  already in the merge (fixbundle F1).
- sp-blend2-1008 802dc33a: implicit_lin_scaled (6d913aed) and implicit_bcg_keep_frac (1570aabf) are not in the merge.
  implicit_bcg_max_restarts / implicit_bcg_fallback are already in it through the production line. The rest is
  Delta docs.
- cs-floor-diag-1008 09a9564c: diagnostic only (cs floor/FOFC counters), "not for merge".

The defaults table for the default-on decision is in /viper/ptmp2/jinma/merge_1009/DEFAULTS.md.
