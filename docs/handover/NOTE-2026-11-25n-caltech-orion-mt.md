# NOTE 2026-11-25n caltech -> orion: MPI / fence timer arms (your 25n; mt files under caltech_results/MT/<arm>/)

## HG (8 ranks, 02-03 + 02-09, HGN's command, binary 79c1cabef 9bebf6a3, libmpit.so built with the hpcx 2.17.1 mpicxx / g++ 13.2.0; no fences added)
- hg mt_base_6: 2026-11-25 02:49:00 rc=0 ended: Terminating on cycle limit wall=106 s CW 0 FATAL 0 mt files 8; user.hst 6 rows vs KT base of the same length: printed-bitwise (6 common); user.hst 6 rows vs base without timer: printed-bitwise (6 common); hydro.hst 6 rows vs KT base of the same length: printed-bitwise (6 common); hydro.hst 6 rows vs base without timer: printed-bitwise (6 common)
- hg mt_base_26: 2026-11-25 02:51:12 rc=0 ended: Terminating on cycle limit wall=131 s CW 0 FATAL 0 mt files 8; user.hst 26 rows vs KT base of the same length: printed-bitwise (26 common); user.hst 26 rows vs base without timer: printed-bitwise (26 common); hydro.hst 26 rows vs KT base of the same length: printed-bitwise (26 common); hydro.hst 26 rows vs base without timer: printed-bitwise (26 common)

## sp stars (4 ranks, 02-16; 25b rsts; MPI / fence timer; 6 and 26 cycles; own keys as production)
- t9l62 mt_base_6: 2026-11-25 02:48:35 rc=0 ended: Terminating on cycle limit wall=81 s CW 0 FATAL 0 mt files 4; user.hst 6 rows vs KT base of the same length: printed-bitwise (6 common); hydro.hst 6 rows vs KT base of the same length: printed-bitwise (6 common)
- t9l62 mt_base_26: 2026-11-25 02:50:52 rc=0 ended: Terminating on cycle limit wall=136 s CW 0 FATAL 0 mt files 4; user.hst 26 rows vs KT base of the same length: printed-bitwise (26 common); hydro.hst 26 rows vs KT base of the same length: printed-bitwise (26 common)
- m35mms mt_base_6: 2026-11-25 02:51:28 rc=0 ended: Terminating on cycle limit wall=35 s CW 0 FATAL 0 mt files 4; user.hst 6 rows vs KT base of the same length: printed-bitwise (6 common); hydro.hst 6 rows vs KT base of the same length: printed-bitwise (6 common)
- m35mms mt_base_26: 2026-11-25 03:01:12 rc=0 ended: Terminating on cycle limit wall=583 s CW 0 FATAL 0 mt files 4; user.hst 26 rows vs KT base of the same length: printed-bitwise (26 common); hydro.hst 26 rows vs KT base of the same length: printed-bitwise (26 common)

### m35mms mt_base_26: anomalous speed, results clean (hst printed-bitwise)
- wall 583 s for 26 cycles, about 20 s per cycle (1.03 s normally; mt_base_6 with the same preload 35 s total, 1 s per cycle). GPU utilisation 0 % on all 4, the node had no other GPU work, load 5.5.
- gdb of two ranks while it was slow: rank 0 in PMPI_Waitall <- CleanHaloFinish <- CleanGhostExchange <- CleanGhostFillT (UCX progress); rank 3 in ioctl <- ibv_reg_mr <- mlx5_reg_mr <- uct_ib_reg_mr <- ucp_memh_register <- ucp_datatype_iter_mem_reg_single (length 4128896).
- mt.0.txt: Waitall totals vet_vgd_ring 192 s, cln_hd_unpack 161 s, cln_hd_pack 63 s, cln_apply_st 24 s (26 cycles). So the time is in registering memory for transfers (the same 4 MB buffers registered again, by the stack), on one node only, and did not happen in any other arm (mt_base_6, KT, ABBA) or in production (1.03 s). No variant tried. It may be the trigger to look for in a slow run of this star.
- Other 25n arms as posted.
