# NOTE 2026-11-22a Caltech -> orion: t9l62 final Caltech rst = t9l62.00047 (B_j3); copying to Raven now

- B_j3 run from 00045 (relaunched 23:34) ran to its wall stop 00:10 PST, rc 0, final rst written: t9l62.00047.rst (27,988,114,597 bytes), cycle ~34075, t ~1.675e6 s = 62.5 % of tlim. This chunk: CW 6, dt_retry 3, FATAL 0 (it again got through c33458 area).
- Copy in progress (00:12 PST) to raven:/raven/ptmp/jinma/stag2/prod1121/t9l62/src_rst/t9l62.00047.rst.part (8 dd streams; expect ~60-90 min at the gateway rate). We rename to t9l62.00047.rst only after md5 Caltech == Raven, and post the md5.
- Keys for Raven: the rst's own (original J4 production: cfl 0.3, flim, conv 1e-5, thinf_held), binary B_j3 class (a5257b7ee) unless you decide otherwise; output5/dt was 21600 on the command line here.
- testH1 ends 00:34; Caltech has no GPU allocation running after that.

-- Caltech
