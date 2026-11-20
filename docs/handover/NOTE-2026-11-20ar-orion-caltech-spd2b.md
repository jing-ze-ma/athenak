# NOTE 2026-11-20ar orion -> caltech (USER): SPD2 gate, P5 step 2 (new build spd2b), medium priority

New build **B_spd2b** = caltech-runs/spd2b, commit 170812c477d6b2b7be30eebcb91e1d40c2e6673f (v20l plus one commit: fewer host syncs; it should be bitwise). Plain cmake line.
Arms: `docs/handover/caltech_runs_1121/arms_spd2b.txt`. Xchecks to put in the result lines:
- SD_rsg hst == SDT_b_1 (B_v20l, 10 cycles), all user columns plus hydro.
- SDR_rsg: rt 3 + 3.
- SDP_rsg (CLEAN_POISON_GHOST=2) == SD_rsg.
- SDT_s_hg == SDT_b_hg.
- Timing: SDT_b/s rounds 1 and 2, as for SPG (together, then swapped). The HG pair is only if cheap.
Priority: after the v20l gate and the ELEAK AB10 arms.

-- orion (for the USER)
