# NOTE 2026-10-21 caltech -> raven: RSG E restarted FRESH on a faster, bitwise-identical build 63a7d28b (USER)

- Branch vetgd-xchg-fast-1021 (Caltech, off 4de8168d; not pushed yet):
  d9b1e07e vet_gd cs band buffers in pinned host memory; 3d624870 opt-in vet_gd_cs_xgroups (default 1 = off);
  63a7d28b VetGdCsCompact builds its lists straight from the packed 2-bit mask (zero words skipped), marks sent as packed words.
- Job 5227122 (2 x 4 H200, prod deck, nlim 50, GDR off, same job): 4de8168d 2.244 s/cycle vs 63a7d28b 1.870 (-16.7 %);
  ic_rad init 155 -> 87 s; list build per rotation 18.0 -> 3.7 s; MPI wait 27.4 -> 8.8 s. All 10 files byte-identical
  (6 bins, hydro.hst, user.hst, 2 rst). GDR on still crashes, now in ImplicitHaloExchange (8 MB RDMA_READ), so GDR stays off.
- USER: start RSG E fresh on it. Chain 5227533 -> 34 -> 35 (2 x 4 H200, 7 d), dir rsg_E_1020/prod63, deck unchanged (a7e734c4),
  binary athena_gpu_he_star_m1_63a7d28b_nofma md5 c0d0e76f. The old chain 5223807 (4de8168d, at ~22 d) runs until 5227533 starts, then is cancelled.
- The UCX reproducibility env (NOTE 6f846a1a) goes into run63.cfg UCX_ENV, read at each link start, once the screens pick one.
- Page https://claude.ai/artifact/LBASg2iK2ArcLY1V2yBCAR now built from the Caltech run (DeltaAI pipeline d316b8a0); it switches to prod63 when that starts.
