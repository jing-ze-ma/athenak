# TASK 2026-10-19 raven -> caltech: 3-D accretor PRODUCTION v11 (stream on), H200 — queue alongside Raven

USER 10-19: the 3-D accretor goes stream-on with v11 and has PRIORITY among the 3-D production runs; "queue it on caltech also".
Raven chain 31274287/88 (2 nodes x 4 A100) is queued with the same setup. USER rule for the accretor (10-18): KEEP BOTH copies
(do not cancel either site without the USER). Never mix sites within one run (different layout -> restart not bitwise).

What changed vs your v3 (v10, stopped per c671ac40): stream_sector_light = neck + donor_irr_scale 0.35 (held light; matches the
analytic donor-light budget; v10 over-irradiated 3.5-8.7x). The onset crash seen on Viper was the MI300A halo bug, not physics:
on A100 the onset to contact runs clean (v11: 16 retries, min dt 1.2e-7, NC 0; run-to-run + restart BITWISE).

1. Build: branch origin/acc3d-held-1018 @ 23e06019 (= v10 a6624af6 + held light), H200, your nofma recipe (as v3). Report md5.
   (A merged production binary prodmerge-1019 may follow; we will tell you — start on 23e06019.)
2. Deck: docs/handover/acc3d-prod-1018/acc3d_prod_v11.athinput.in (replace @FILES@ with your copy of this directory; column + both
   tables unchanged from v3: md5 5b813886 / d81b7a46 / 99d668f6). Layout as v3 (2 H200), tlim 1e9 in the deck.
3. USER rule before production (10-19): on the H200 build, a short determinism + restart check: from t=0 with
   problem/stream=true and a short static-free start is NOT the production path, so do it like this: run the static phase to
   t ~0.002 (C), rerun identically (C2), and restart C from a mid rst (R) -> C vs C2 and R vs C bitwise (hst every cycle + final bins).
   If DIFF: report, do not start production.
4. Production: your acc3d chain (acc3d_chain.sh: static phase to t 0.05 with problem/stream=false, then the deck as is, stream on),
   FRESH from t=0, 12 h links, as many links as v3 had queued. STOP only on FATAL / DIVERGED / summary NON-CONVERGED=[1-9] / NaN / rc.
5. Report: md5, check results, link ids + start estimate; after link 1 starts: s/cycle, dt, memory/GPU, Picard summary.
