# TASK 2026-11-09 viper -> Orion: CPU batch (USER: "anything that can be run on cpu can go to orion")

Run these self-contained TASKs in this order (all CPU, MPI build lines inside each). Answer each with its own NOTE
as the TASK says. Independent of each other; run in parallel if Orion has the cores.

1. `docs/handover/TASK-2026-11-09-orion-hgpos.md` - HG positivity: time2_adm_eflr + implicit_pos_transport_t gates
   (fork/hgpos-1109b b798f1a30 vs prod-1112-cand f34ada236). NOTE: NOTE-2026-11-09-orion-hgpos.md.
2. `docs/handover/TASK-2026-11-09-orion-r136fast.md` - R136 floor-gas cap + re-block identity gates
   (fork/fastfloor-1109 138e6ecb5). NOTE: NOTE-2026-11-09-orion-viper-r136fast-gates.md.
3. `docs/handover/TASK-2026-11-09-orion-raycheap.md` - physical-clock ray-average gates + 2nd-order convergence check
   (fork/raycheap2-1109 a39c9ec49 vs prod-1112-cand 3330bc0dd). NOTE: NOTE-2026-11-xx-orion-raycheap.md.

Not included: TASK-2026-11-08-orion-hg-mesa-lin.md (USER declined on Orion, see NOTE-2026-11-08-orion-viper-hg-mesa-declined.md).
4. `docs/handover/TASK-2026-11-09-orion-lintol.md` - R136 small (831x8x8) Picard stop-test accuracy: find the cheapest
   stop setting reaching ~1.5e-8 vs the reference (deck/data in docs/handover/lintol-1109/; fork/prod-1112-cand).
   NOTE as the TASK says.
