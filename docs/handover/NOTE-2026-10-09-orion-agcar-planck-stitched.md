# NOTE from Orion: stitched AG Car Planck table ready (follows NOTE-2026-10-09-orion-lowrho-planck-done)

The user approved the anchored stitch proposed in the -done note. Orion built it.

- Commit 4dc8f1d9 on `agcar-opac-orion-1009`.
- Table: `docs/handover/agcar-opac-orion-1009/stitched/planck_ext2orion_gs98_x0.36_z0.02.txt`, md5 493f9d01.
- Format: the same as ext2 (dc482452): 91370 lines, the same grid line, read with he_planck_table
  (HsReadOpacityTable) unchanged. Two header lines (3 and 9) are reworded to describe the stitch.
- Rosseland: use ext2 `rosseland_ext2_gs98_x0.36_z0.02.txt` (b28e97c5) as is.

## What changed
- 4470 cells, all at log T 3.45-4.30 below ext2's data floor rho_f(T). rho_f(T) comes from the build_ext.py fergR
  logic and is listed in `stitched/floors.txt`:
  - log T <= 4.0: the log R = -8 line.
  - log T >= 4.2: the TOPS floor, log rho -15.
  - log T 4.025-4.175 (the blend): the Ferguson floor.
- The other 86887 cells are text-identical to ext2.
- There, log kP = ext2(T, rho_f) + [orion(T, rho) - orion(T, rho_f)].
- The correction is weighted 1/3 and 2/3 at the two T nodes inside each edge, and 0 outside 3.45-4.30.
- The correction ranges from -1.62 to +1.11 dex at log rho -21.
- Continuity: the jump across rho_f is at most 0.039 dex. The largest jump between adjacent T rows is 0.68 dex at
  log rho -21, against 0.57 in ext2.
- Details: `stitched/STITCH.md` and `stitched/checks.txt`. Plots: `stitched/stitch_delta_map.png` and
  `stitched/stitch_cuts.png`.

The Planck mean above 11 kK is anchored to TOPS/blend at the floor, so the orion normalisation deficit does not
enter.
