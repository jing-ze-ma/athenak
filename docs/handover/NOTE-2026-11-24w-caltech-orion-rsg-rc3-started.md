# NOTE 2026-11-24w Caltech -> orion: RSG running on rc3 from rsg.01707 (3 nodes)

- rst rsg.01707.rst: transferred over ssh+zstd in 38 min. md5 c2da5871810bab9f9af77809911794de equals viper's. The 3 table files were copied too (md5s equal).
- Run: rsg_rc3_1124/, holder 7333072, nodes 01-17 + 02-03 + 02-09, 12 ranks, step .89.
  - Binary athena_ia123rc3_af836f5a3 (09f50abe) with IA123_rsg_restart.args (already has fofc + fofc_retest + mg_levels 4).
  - Table paths repointed: problem/he_ic_file, he_opac_table and he_planck_table point to the Caltech copies.
- The first try (step .88) died before cycle 0: "ib_mlx5 Local protection error" (RDMA_WRITE) in UCX on the first cross-node message.
  - Fixed with our old 10-21 RSG multi-node env: UCX_TLS=^xpmem UCX_IB_GPU_DIRECT_RDMA=no UCX_RNDV_SCHEME=put_ppln (put_ppln was run-to-run bitwise on 2x4 H200).
  - Worth adding to the RSG args notes / release docs.
- First cycles: 0 FATAL, 0 CW, 0 created inf. Cycle 771891 -> 772000 (t 7.2241e8 -> 7.2251e8 s), dt 1000.0 (cap), ~1.54 s/cyc.
  - Last hst: t 7.225340e+08 L_top/L_in 0.92382 L_bot/L_in 1.04518.
  - fofc_retest retested 6.0e8 cells.
- Wall to the holder end is ~51 h.

-- Caltech
