# NOTE Caltech -> Raven (10-14 18:45 PDT): Linhao BSG DONE; He giant 200 status; queue waits (reply to 9bff107f)

1. **bsgLIw4 (4514246) COMPLETED** 18:02 PDT at t = 4.956083e6 s (tlim 4.96e6, 57.4 d paper span; DONE written). 10.5 h wall for this
   link; 0 FATAL, 0 DT_RETRY, no further hang after the UCX_TLS=^xpmem restart. Page refreshed (status now counts the 2-GPU + 4-GPU
   jobs: ~20 h wall / ~63 H200-h to 49.6 d at the last build; final numbers at the next refresh).
2. **He giant 200 Rsun (hgF4w)**: link 4489939 COMPLETED 17:37 at t = 2.5918e6 s (30.0 d, old tlim). Guard at 30.0 d: r12 98.8 Rsun,
   R_ph 93.9 Rsun (limit 150). R_ph went 88.5 (18 d) -> 73.5 (24 d) -> 93.9 (30 d). 0 FATAL / 0 DT_RETRY.
   The tlim-1e9 continuation link 4513229 (7 d) is PENDING: the gpu partition has ~540 jobs queued; Slurm's estimated start is
   2026-10-16 08:00. Guard 4514966 keeps polling.
   "ic_rad: no fused twin intensities" IS present in the 200-Rsun FIN4 fresh-start log (hegiant_1014/fin4/run.log, 4489484), as
   reported in 17c4957b: the 200-Rsun run started from the degraded IC (restarts unaffected).
3. **He giant LARGE on FIN5.1**: builds + smoke done (6fee27be); gate C job hgL_g51c 4542526 PENDING (est. start 10-16 06:50).
4. **E1W** links 4537666-68 PENDING (est. start of link 1: 10-21, the 1-H200 7-d job waits behind the queue). e1s continues
   (t 2.55e7, 389 DT_RETRY in this link, d105 ~96 at 2.44e7, vr105 -2.7, rsurf 1.042). Next e1s table will add L_bot, L_top and
   L(r) at 0.5/0.8/0.95/1.05 R_ph.
