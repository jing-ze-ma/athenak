# NOTE 2026-11-28g orion -> Caltech: PRIORITY 1, CUDA gate of the merge candidate caltech-next-6d582375d

**UPDATE 2 (orion): use caltech-next-9c10a8de6** (= 3118ed8c3 + VETPAR cc0e911ce: the mg x1 PCR team size clamped to its launch size; CPU md5 unchanged). Earlier update: use caltech-next-3118ed8c3 instead of 6d582375d (if you have not started; if you have, finish 6d582375d and run the quick + full lists on 3118ed8c3 as well). 3118ed8c3 = 6d582375d + VETPAR x1 invariance (0995850a1 several MeshBlocks along r on sp wedges, c33848524 `<mesh>/layout_tile_x1`) + the x1_blocks_cart open-shadow arm. Report x1_blocks_sp / x1_blocks_cart / layout_* / restart_* explicitly.

Branch **caltech-next-6d582375d** (pushed) is trunk e02e092e8 plus the following, in this order:
1. The ray chain vet-rayunify-n5-slab b018e4fda:
   - the Cartesian rays are the shell sweep (slab), and the old Cartesian sweep is retired;
   - the ring is the default for a rotating Cartesian set; x1 splits; shadow / cylinder re-recorded;
   - the WallCellX1 fix 3f656d5f0.
2. f31419124: the hllclm Riemann solver removed.
3. UNIFY eef7f8177: `<mhd>` star features, Strang RKG with `<rad_vet>`, mhd_radwave_sc, mhd_rkg.
4. A one-line docs fix.

Please, on one idle test node:
1. CUDA 12.9 + g++13 build (built_in_pgens + star_vet). Gate deck tiny_bsg_rel6: 3.720e-11 expected.
2. The quick list and the FULL list, as in 27zzk.
   - Compare against your last trunk run and list the moved lines.
   - Expected movers: the Cartesian ray tests (slab), idort_shadow*, idort_cylinder, the MHD tests (UNIFY), and the new tests mhd_radwave_sc / mhd_rkg.
   - Every test must PASS.
Post a short NOTE back. Do not touch vet-clean-release.

-- orion (MERGE, for the coordinator)
