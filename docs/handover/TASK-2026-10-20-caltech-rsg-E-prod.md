# TASK 2026-10-20 raven -> caltech: queue the FRESH global RSG design E on 8 H200 (third race copy; USER "queue")

The fresh global RSG (design E) is also queued on Raven (6x4 A100) and Viper (12x2 MI300A); first to START wins, the others are
cancelled only after the winner has started (Raven tells you). Your copy is the long-shot backup (est. ~96 h); queue it now.

1. Build: origin/prodmerge-1019b @ 563ab474 (= fin51-r-fix aaa01b3f + time2_reject=retry + relaxed_esc held-gas fix), he_star_m1,
   H200, your usual nofma recipe (as for rsg g25). Report md5.
2. Deck: docs/handover/rsg-E-1020/rsg_E_start1020r.athinput.in (Raven original md5 ae8d2e36). Replace @FILES@ with your copy of
   docs/handover (the deep IC column docs/handover/rsg-E-1020/ic_rsg_g25_deep.txt md5 2d7de3c8; Rosseland/Planck tables in
   docs/handover/rsg-g25-prod/, md5 fa0f4f42 / a9a8d6d3, unchanged). Design E: nx1 548 (geometric head), merged lateral cells
   f 8/4/2 below 0.1 R_ph, cs N=128 full sphere, 24 MeshBlocks 548x64x64 -> 8 H200 (3 blocks/GPU; check memory: ~25-30 GB/block).
   Key method: vet_source = relaxed_esc, implicit_am_cs = true, time2_reject = retry (+ dt_retry_fac_t2 0.9, grow 1.1, ngood 20),
   cs_seam_flux = shared (default). Fresh start, tlim 1e9 (walltime-limited links).
3. Quick check before the chain (USER rule): reduced N=32 on 1 node: C, C2 identical + R (restart) == C (hst + bins), as for rsg g25.
4. Chain: your RSG chain script (as for rsg_g25_prod), FRESH from t=0, 7-day links or your max, afterok on the check.
   STOP only on FATAL / DIVERGED / summary NON-CONVERGED=[1-9] / NaN / rc. Report link ids + start estimate.
5. When Raven says another copy has started: cancel yours (keep anything already written).
