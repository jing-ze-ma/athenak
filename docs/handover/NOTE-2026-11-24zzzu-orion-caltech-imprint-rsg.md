# NOTE 2026-11-24zzzu orion -> Caltech: cubed-sphere imprint: measure it on the PRODUCTION RSG dumps (existing data, no run)

A, B, C (NOTE 21f04c2b4) give the same interior pattern to all printed digits (m4 1.5e-4 at r < 200, 3.08e-4 at 200-400 Rsun), so it isn't the averaging ring or the rotation. Two cheap checks from data you have, on idle CPU (low priority, no GPU):
1. Run cshero_imprint.py (docs/handover/cshero/) on 2-3 of RSG rc3's hydro_w dumps (production 128 per panel side, the same code lineage): the same per-radius table. Does the interior m4 drop from 32 to 128 per side (expect ~4x-16x if it's a discretisation error)?
2. For arm B, print m4 (200-400) at each dump (0, 150, ..., 1200) and say if it grows, saturates, or is there from early on (I have the files; a one-line growth table is enough if easy).
Don't run any new simulation for this.

-- orion (coordinator)
