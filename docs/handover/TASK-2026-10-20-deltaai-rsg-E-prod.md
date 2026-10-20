# TASK 2026-10-20 raven -> deltaai: queue the FRESH global RSG design E on 2 ghx4 nodes (8 GH200) as a race copy (USER "put it on delta ai")

Race: the same run is queued on Raven (est ~19 h), Viper and Caltech; first to START wins, Raven tells the others to cancel
(cancel only after the winner has started). Your budget (867 SU, ~108 node-h) covers ~2 d on 2 nodes: if you win, run until ~100 SU
remain, then stop cleanly at an rst and say so (we continue elsewhere from your rst).

1. Build: the NEWEST origin/fin51-r-fix at build time (now 0fbfa31f; 82275582 = + retry BE fallback, torque registers, hr_q opt-in,
   is being gated and will be pushed soon — if it is on fin51-r-fix when you build, take it). he_star_m1, GH200 (Kokkos_ARCH_HOPPER90),
   your nofma recipe (-ffp-contract=off), MPICH_GPU_SUPPORT_ENABLED=1. Report sha + md5.
2. Deck: docs/handover/rsg-E-1020/rsg_E_start1020r.athinput.in (as Caltech: replace only the @FILES@ lines; IC ic_rsg_g25_deep.txt
   2d7de3c8, tables in docs/handover/rsg-g25-prod/ fa0f4f42 / a9a8d6d3). 24 MeshBlocks 548x64x64 -> 8 GH200, 3 blocks/GPU.
   Fresh start, tlim 1e9.
3. No pre-start check job (USER: skip unneeded tests). Put the run + rst on /work/hdd/bivj (nvme is nearly full); prune rst to 3.
4. Chain: 48 h links on ghx4 (-A bivj-dtai-gh, 2 nodes, 4 tasks + 4 GPUs per node), afterok on the previous link, at most as many
   links as the budget allows. STOP only on FATAL / DIVERGED / summary NON-CONVERGED=[1-9] / NaN / rc.
5. Do NOT touch the ghx4-interactive slot (it belongs to the BSG chain).
Reply in a NOTE deltaai->raven: build sha/md5, job ids, start estimate, and the ghx4-interactive limits (max walltime, nodes/GPUs per
job) for reference.
