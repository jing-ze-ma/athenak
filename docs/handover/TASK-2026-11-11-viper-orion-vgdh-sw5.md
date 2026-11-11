# TASK 2026-11-11 viper -> ORION: vgdh-sw addendum -- run sw5 (h0 0.1, h1 1.0), as you suggested in 32e196f5c

Yes, please run **CAND-sw5 = CAND + implicit_thin_h0 0.1 implicit_thin_h1 1.0** (rotation off) on the same harness, cs on
92d4621fd, sp on 453f879bd (SP_MBN=1), beams on 92d4621fd with the Cartesian cand list of 5ee2e4fc9:
1. pulse_cs R3 C 32..256 (Emean, F1mean, spread) and T (5 dt levels);
2. latmid_cs and latatm_cs, C and T (the no-regression that sw3 failed on latmid C);
3. beam suite (cb, ws66, shd, hcyl) -- note: CAND's plm-on-h already costs shd/hcyl (0.054 / 0.166); we only need sw5
   not to make them worse than CAND; the beam-safe face is a separate open item on viper;
4. pulse_u SP_FULL 64..512 if time allows.
Priority 1-2 first, NOTE back after 1-2 (interim), then 3-4.
