# NOTE Caltech -> Raven coordinator: status 10-13 11:50 PDT (please read)

1. **AG Car E1: READY, waiting only for your GO.** Smoke 4454178 GATE PASS (C vs L 30+30 and C vs Coff all BITWISE, 8/8 dumps,
   11/11 hst rows; dt 727.27 every cycle; Picard 2.39/6; NC 0; 0.60 s/cycle on 2 H200 -> about 5.6 h, 1 link). Prune test 6/6.
   Links 4452279 (HELD) -> 4452280 -> 4452281. **Push docs/handover/GO-e1-caltech.txt and link 1 is released within ~15 min**
   (start estimate was ~14:00 PDT). Details: NOTE-2026-10-13-caltech-e1.md.
2. **Two bugs in the e1caltech-1013 bundle (worked around here; please fix for other sites):**
   - smoke_e1.sh writes `${LAUNCH:+LAUNCH="$LAUNCH"}` through an unquoted heredoc, so a multi-word LAUNCH loses its quotes
     (first smoke 4452278 rc 127, nothing ran). Workaround: one-word LAUNCH script. Fix: escape the quotes, or say "LAUNCH must be one word".
   - cmplast.py / chaincmp.py -> bin_convert.py import h5py at top level; not installed by default here (pip --user fixed it).
3. **He giant FIN3: RUNNING** (link 1 4454171 since 11:32 PDT, chain 4454171-76). Smoke 4452306 clean, restart BITWISE, 85.7 GB/H200.
   **Speed 2.74 s/cycle on 2 H200 x 2 blocks (not ~1.4)** -> 30 d needs ~9 x 12 h links, not 5-6. 6 are queued; more only with
   the user's OK after link 1 measures the real rate. If another site is much faster, say which copy should carry it.
   Details: NOTE-2026-10-13-caltech-hegiant-fin3.md.
4. **BSG fin3:** link 1 4434961 near its 12 h limit (t 1.35e6 of 4.96e6 at 06:00; s/cycle drifted 1.2 -> ~2.3), links 2-3 queued;
   may need a 4th link. Per-link report goes to NOTE-2026-10-13-caltech-bsg-fin3.md when link 1 ends.
5. **od mbband + open-ot ENTHBE:** received, in progress (mbband first, as asked). QUEUED notes follow under the TASK names.
