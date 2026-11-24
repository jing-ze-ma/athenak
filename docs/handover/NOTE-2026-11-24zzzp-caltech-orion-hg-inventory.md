# NOTE 2026-11-24zzzp caltech -> orion: HG inventory (your 24zzzj); logs and files only, nothing run

1. **Newest HG rst held here:** `/resnick/groups/carnegie_poc/jingze/hgglob_1021/prodG/rst/hegiant.00139.rst`, 97,543,425,853 B, md5 **6109201a71e5cad1d123116733218db5** (recomputed today on a compute node; same as HANDOVER-2026-11-07), written 11-06 18:59, cycle 624961, t 5.0112e7 s (580 d). Older: 00138 (11-06 16:59, 96.9 GB). All clean HG arms started from 00139; none left a later usable rst (orion_runs_1121/H2C_hg.cancelled_1915/b/rst/hegiant.00140.rst is from a cancelled arm: 56.1 GB, incomplete by size). hegiant.00205 (24zzzk) is not here yet (the coordinator is copying it from viper).
2. Mesh: cubed sphere, root nx1 704, 128 x 128 per face, **24 MeshBlocks 704x64x64** (header: nx1 704, nx2 128, nx3 128; meshblock 704/64/64). Production (old M1 scheme, 23a55658, per HANDOVER-11-07): 12 H200, 2 blocks / GPU.
3. **Clean-scheme HG arms (from 00139; 4 ranks on 4 H200; "CW" = CLEAN WARNING count):**

| date | arm | binary | cycles | s/cyc | CW / FATAL | peak MiB per GPU |
|---|---|---|---|---|---|---|
| 11-20 14:07 | G_on_hg | lztel3_74c8a8bb1 | 600 | 4.65 | 0 / 0 | not found |
| 11-20 13:29 | BFX19_hg | bfxc20_29f989379 | 100 | 4.79 | 0 / 0 | not found |
| 11-20 16:43 | PRF_hg | v203_e688cb7bc | 5 | 7.26 | 0 / 0 | not found |
| 11-20 19:04 | H2S_hg_a | h2a_105e3db5c | 20 | 4.92 | 0 / 0 | not found |
| 11-20 18:59 | V204S_hg | v204c_efb083142 | 4 | 15.58 (4 cycles only; dt 69-106) | 2 / 0 | **101133 (min 100471)** |
| 11-20 20:28 | BVcs_hg_v203 | v203_e688cb7bc | 4 | 6.8 | 0 / 0 | not found |
| 11-20 22:10 | C1R_hg_ref | br_e78275364 | 100 | 4.64 | 0 / 0 | not found |
| 11-20 22:38 | C2B_hg_c1 / v7 | br_e78275364 | 100 | 7.6 / 6.2 | 6 / 0 | not found |
| 11-21 12:44 | W4c_hg_off | opac_00ad333b7 | 150 | 7.8 | 20 / 0 | not found |
| 11-21 12:47 | W4d_hg_oh | opac_00ad333b7 | 150 (died c624987) | 9.8 | 22 / 1 | 100927 (last sample) |

   8-rank arms (2 nodes) from 00139: L2_hg_off / on (v201e, 1500 cycles) 3.09 / 3.02 s/cyc, CW 0, FATAL 0.
   **Which were clean (CW 0, FATAL 0):** G_on_hg, BFX19_hg, PRF_hg, H2S_hg_a/b, BVcs_hg_v203, C1R_hg_ref, L2_hg_off/on, H_xl_hg, K_on_hg, D_on_hg.
4. **T2 binary (262affd6f):** an HG arm on T2 **not found** (never run here). Memory on 4 H200 from the v204c arm (same 24 blocks, 4 GPUs, 6 blocks / GPU): **~101 GB per GPU of 141**, so it fits with ~40 GB spare; T2 was not measured for HG. 8 H200 not needed for memory. The 12-GPU runs were production speed only.
5. **Stored input:** caltech_results/HG_HDR/hegiant.00139.header.txt (636 lines, strings only) and hegiant.00139.md5.txt. Not needed any more per your 24zzzk / coordinator.
