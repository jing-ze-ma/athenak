# NOTE 2026-11-24zzzg orion -> Caltech: fixes to the CSHERO NOTE 24zzn (arm C deck was broken)

1. **Arm C:** star3d_cs_C_fixed.athinput was wrong (its own Orion run FATALs at init: clean_rot_ring = live with rotation off). Use **docs/handover/cshero/star3d_cs_C_fixed2.athinput** instead. It hasn't been run either: run `time/nlim=0` first; if it FATALs, skip C and post the message (don't repair it yourself). C differs from B in rotation AND the ring / cross-fade / EMA keys.
2. **Text corrections to 24zzn:** that the hero ran the RELEASE template is inferred, not established (no Raven log or dumps found); the K4-window cause is a hypothesis; don't read absolute sym_rms values as imprint: compare arms.
3. Run order unchanged (A, B first; C2 last, only if B keeps the X or as a control). Before starting, run `time/nlim=0` on A and B too: their decks haven't run on your binary.

-- orion (coordinator)
