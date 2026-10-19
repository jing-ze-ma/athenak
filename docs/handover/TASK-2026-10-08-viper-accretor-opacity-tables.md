# TASK 2026-10-08 for viper: X 0.70, Z 0.02 opacity tables for the accretor RHD run (Caltech)

From: Caltech (accretor owner). User 10-08: "viper has new opacity tables, ask viper".

## Need
Rosseland AND Planck mean opacity tables for the Plaskett gainer composition, **X 0.70, Z 0.02** (GS98 mix if you have it;
say which mix and source you used), for the accretor RHD port (branch accretor-rhd-1008, ry_per_accretor + rad_m1, tables
read like he_star_m1's he_opac_table / he_planck_table).
- Range wanted: log T 3.5 (or your table floor) to ~7.6, log rho -18 to 0 (g/cc). If your new tables stop lower in T
  (the older TOPS ones stop at log T 7.065), say so: only the 300 km/s hot ambient (log T 6.83, rho ~5e-16, absorption
  masked) is above ~6.5, so log T 7.0 is acceptable.
- Format: the same text format as the He giant / AG Car tables (rosseland_tops_*.txt / planck_tops_*.txt, as read by
  he_star_m1 / rad_m1 table opacity); state the grid (dlog T, dlog rho) and any edge handling.
- Please push them with md5s to `docs/handover/accretor-tables-1008/` on **accretor-rhd-1008** (fetch, add, push; never
  force) plus a short NOTE (source, mix, range, md5; if they are your "new tables", what changed vs the TOPS ones).

## Why
S3a (radiation on the static Roche envelope) passes all gates with the dev pair TOPS X 0.7 Z 0.008 (bsg_1001 bundle),
but the run should use the real Z 0.02 (stronger Fe bump). Caltech will rebuild the IC column
(docs/dev/accretor_rhd/make_ic_accretor_column.py) and rerun S3a with your tables. Status: docs/dev/accretor_rhd_design.md
section 9 on accretor-rhd-1008.
