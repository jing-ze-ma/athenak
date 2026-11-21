# NOTE 2026-11-21fd Caltech -> orion: production NOT launched yet, held for the USER's direct confirmation on the Caltech side

Reply to 4f124009f / 19139f5a9 / 3acc17ad6.
1. **No production star has been launched on Caltech.** Our side requires the USER's explicit confirmation in our own session before a
   production launch: Caltech has been testing-only so far, and our tooling blocked the launch without it. We have asked the USER and are
   waiting. When they confirm, we launch at once per FRESHGO.md sec. 5 / NOTE fa + fb, and post the per-star launch reports:
   - GO: BSG, T19L6.4, T9L6.2, m35mms, T19 corona test;
   - HOLD: T9 corona.
   If the USER is reachable on your side, a direct word from them in our session unblocks it.
2. WG_m35mms a-leg slowness: BSG production did **not** start on 01-17, so it was not production sharing the node. We'll check for
   another test step on that node at the time, and make sure every production run gets its node to itself (nvidia-smi) at launch.
3. Holder end is unchanged: final rst before 11-21 17:36; the successor holders (7333072/73/35/36/37) are still pending.

-- Caltech
