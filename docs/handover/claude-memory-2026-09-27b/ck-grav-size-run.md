---
name: ck-grav-size-run
description: ck_grav_size -- the OVERSIZED sizing run; ANSWERED 2026-08-26 at 21 rotations, r99 = 1.794e10
metadata:
  type: project
---

`/viper/u2/jinma/ATHENAK/bench/ck_grav_size`, set up 2026-08-26. Binary built fresh from
HEAD **c37ebe75** (first commit with `problem/grav_point_mass`); input is `ck_limb`'s
verbatim except the four rows below. Full rationale and the measurement recipe are in that
directory's `NOTES.md`.

| parameter | ck_limb | here |
|---|---|---|
| `problem/grav_point_mass` | absent (false) | **true** |
| `mesh/x1max` | 1.433e10 | **1.90e10** (deliberately oversized) |
| `mesh/nx1`, `meshblock/nx1` | 88 | **172** (holds dr = 5.558e7) |
| `time/tlim` | 8.64e7 | **4.32e7** = 500 days = 141.6 rotations (raised 2026-08-26 from 1.22e7) |

**This is NOT the science campaign.** It answers one question: where does the 1e-6 bar
isobar settle at lon = +-90, |lat| < 30, after the atmosphere relaxes? Then trim x1max to
that + ~3 H, set nx1 to hold dr, and launch production. See [[grav-point-mass-flag]] for
why both existing estimates (1.64e10 and 1.87-1.92e10) are untrustworthy -- both are
near-IC, and [[dhj-floors-for-1e-6-bar]] records the same trap.

Measure on the ISOBAR, never a constant-r shell ([[dhj-isobar-vs-shell]]), and recover p/T
with `<hydro>/eos_table_dump` ([[eos-table-dump]]), never a naive root find
([[eos-inversion-nan-trap]]).

Cost is roughly neutral against [[ck-limb-run]] despite 1.95x the cells, because dt is
2.3x larger with the point mass. Resubmit with `sbatch resubmit_viper.sh` per 24 h timeout.

**Job 11043450 on apu1, submitted 2026-08-26** (queued behind ck_limb 11022482).
Smoke-tested first on apudev (11043426, 300 cycles): runs clean at nx1 = 172 -- no LDS
problem -- at **dt = 17.4-19.2 s and 24.4 cycles/s, so ~12 min per rotation and roughly
8 h for all 40**. dt confirms the 2.3x point-mass gain (ck_limb sits near 10 s on half
the cells).

**Rate, measured at 2.4 rotations: 292.6 simulated s per wall s** (23.8 cycles/s x dt
12.3 s), so 500 days costs **~41 h = two 24 h slots**. dt settled within the first
rotation (19.2 -> 12.3) and is then FLAT; ck_limb behaves the same way (flat +-1 s from
rotation 0.4 to 47.8), so there is no slow decay to price in.

**A RESTART TAKES tlim FROM THE RESTART FILE'S EMBEDDED INPUT, not from the athinput.**
Job 11043450 was launched with the old tlim = 1.22e7 and will therefore stop at 141 days
even though the input file now says 4.32e7; `resubmit_viper.sh` passes
`time/tlim=4.32e7` on the command line to get past it, and job **11046058** is queued
with `--dependency=afterany:11043450` so the continuation starts the moment it exits.
One more slot after that.


## ANSWERED 2026-08-26 at 21 rotations -- the run has done its job

**out% = 0.00 at every one of the 22 dumps**, IC included, so oversizing worked and this is
a measurement. r99(1e-6 bar) over the terminator columns = **1.794e10 cm**; H at the isobar
**1.24e9 = 22 cells**; 0.00 % of isobar cells at either floor.

**H and T converged by rotation 4. Only r99 still moves, and it moves INWARD** (-1.09e7
cm/rot, decelerating), so the buffer at a fixed x1max can only grow -- which is why the grid
could be chosen without waiting for the remaining 120 rotations. An exponential fit says
r_inf = 1.795e10, tau = 4.5 rot, but its residual is 0.7 cell and it under-predicts the late
drift: quote the sign, not the fit.

Production run set up as [[ck-grav-prod-run]] at x1max 2.0556e10 / nx1 200. Method and the
inversion validations: [[isobar-buffer-calibration]].

**This run can be cancelled** (it and its queued continuation 11046058) if apu1 slots are
wanted for production; it is worth continuing only to confirm the inward drift never
reverses.
