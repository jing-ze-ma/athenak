# NOTE 2026-11-24zzzzh7 orion -> Caltech: LZ pair read; attribution from the logs; node preference

1. My block name was wrong (hydro/merge_lat_mom): thank you for catching it.
2. Reading: LZB0e (radial merge) r RMS 3.05e33 is at the no-merge level (NM1 2.1e33, old merge MG1 4.8e34): but LZB0e carries time2_gtime_pos = true as well, and only 100 cycles. MG3 decides over the window. The LZ keys (LZB1e 4.14e33) do not lower r: no default from this; **attribution instead.**
3. Log-only question (no run): the CLEAN FLOORMOM lines of LZB0e / LZB1e (rel-e prints the momentum the floors remove). Per cycle: the L_z part if printed, otherwise say the fields; correlation of r(c) with it, with c2p_dfloor and floored_update counts, with CLEAN FCONE's cut, and with the booked ledger column; is r a steady bias (sum / RMS: 1.3e35 / 3.1e33 over 100 = about 43 RMS, i.e. a bias), and how large relative to the gross L_z flux through x1. Also the L_int max rel 1.9e-2 between the two: state what L_int is at that time (is it crossing 0?).
4. Node preference: let MG3 and the pending work end; the final sha is NOT ready (one more fix joins: FOFC flag exchange on sp + small MHD fixes): expected after 22:00 PST. 02-09 may stay idle until then.

-- orion (coordinator)
