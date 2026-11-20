# NOTE 2026-11-20ac orion -> caltech (USER): LZTEL, cs L_z fix; 5 arms (new build lztel)

New build **B_lztel** = caltech-runs/lztel, commit 5f18abdc7 (v20j plus a new key, rad_m1/clean_cs_lz_tel, default false).
Arms: `docs/handover/caltech_runs_1121/arms_lztel.txt`, all with env LZCS=1:
- T1_rsg: key true, 300 cycles.
- T0_rsg: key false, 300 cycles (the control).
- T2_hg: key true, 50 cycles.
- T1rt_rsg: rt 5 + 5. T1p_rsg: poison, 10 cycles; it must equal T1rt_rsg/a.
Priority: after the queued CONV2/ELEAK short arms, ahead of batch 2/3 remainders. Run T1_rsg and T0_rsg together.
- The key normally sits in the deck. On a restart run it should be taken from the command line ("added from the command line").
  If instead it FATALs with "not found", send the line and skip.
- Send back: user.hst, run.log.gz (the 'CLEAN LZ' and 'LZCS cycle' lines), s/cycle, and any CW / FATAL.

-- orion (for the USER)
