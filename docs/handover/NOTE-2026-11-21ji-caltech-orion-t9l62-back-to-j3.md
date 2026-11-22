# NOTE 2026-11-21ji Caltech -> orion: USER — t9l62 back on the OLD joint3 build from t9l62.00031

USER (after NOTE jh, the 5-10x grid noise on the new builds): "stop it now and evolve it with an older joint3 build from 31".
- B_p1 run on testH1 stopped cleanly 19:3x PST via the STOP file (final rst 00040 written; rsts 00037-00040 of the new-build segments moved to t9l62/rst_p1_seg/, kept).
- Relaunched ~19:35 on testH1 (01-11) from t9l62.00031 (t 15.61 d, the last rst of the B_j3 production, before any new-build key) with athena_j3_a5257b7ee (md5 f71a22b6) and NO command-line keys: the restart file's original production keys (cfl 0.3, flim, conv_max 1e-5, no fwb/fcone/harmonic/opac_update/gas switch). Wall -t to 00:20 PST (B_j3 has no graceful stop, so the last periodic rst is what survives the 00:34 holder end).
- So the Raven hand-off (lf/jf/jg) now means: the newest B_j3 rst after tonight, not 00037. Do not start t9l62 on Raven.

-- Caltech
