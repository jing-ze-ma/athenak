# NOTE 2026-10-20 raven -> caltech: newest production binary = fin51-r-fix aaa01b3f (USER newest-binary rule) — switch at next link boundary

fin51-r-fix moved ab758e94 -> 5864f189 -> aaa01b3f (pushed). Contents since ab758e94 (all gated: strict gbit default path bitwise on
A100 + MI300A + CPU for bsg/e1/hg/hg145; det/restart/retry bitwise incl. 24-rank multi-node):
- multi-node determinism fix (Kokkos::fence before MPI_Irecv into device buffers): 12+-rank runs were NOT run-to-run bitwise before;
- cs shared seam flux (cube-edge instability fix), am_transport for all variants + new FRESH-run defaults (recipe), implicit_am_cs fix,
  accretor wall/tracer/history fixes, cs-s7 + merged-cell options, relaxed_esc (opt-in).
Restarts whose rst lacks the new keys keep the OLD defaults -> running non-cs chains continue bitwise (bin/rst headers gain one
parameter line). cs chains: the seam default changed (legacy 'positive' reads as 'shared') -> a method change for cs restarts.
ACTION: (1) list every Caltech production chain you have (name, binary SHA/md5, deck, next link time) in a NOTE;
(2) for each non-cs chain, rebuild aaa01b3f (he_star_m1 nofma as usual) and switch BIN at the next link boundary; report md5;
(3) cs chains: don't switch without asking (method change); (4) acc3d trc: keep 13206541 (identical results; it's in aaa01b3f's history).
The fresh RSG E copy (TASK 17e12e87) uses prodmerge-1019b 563ab474 = aaa01b3f + opt-in retry mode, as on Raven/Viper.
