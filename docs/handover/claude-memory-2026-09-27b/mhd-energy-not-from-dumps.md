---
name: mhd-energy-not-from-dumps
description: Magnetic energy CANNOT be computed from a bin dump -- bcc is the mean of the face fields while the code integrates the mean of the face energies; they differ by 8.8x in B_r on sp_dhj_ctl
metadata:
  type: reference
---

**2026-09-02.** `history.cpp:341-343` integrates magnetic energy as

    vol*0.25*(SQR(bx1f(m,k,j,i+1)) + SQR(bx1f(m,k,j,i)))

-- the mean of the two face ENERGIES. The binary dump carries `bcc`, which is
`0.5*(bx1f(i) + bx1f(i+1))` -- the mean of the FIELDS. mean(B^2) >= (mean B)^2, so a
dump-based magnetic energy is always LOW, and by a factor that depends on how much the
field varies across a cell.

Measured on `sp_dhj_ctl` at 100 rotations, my-ME / history-ME per component:

    ME1 (radial)  0.113      <- B_r reverses at the GRID SCALE
    ME2           0.529
    ME3           0.680

**So: quote magnetic energy from the history file, and use `bcc` only to show WHERE the
field is.** The 0.113 is itself worth knowing -- it says the radial field is largely
cancelling between adjacent faces on this run.

## The gate that caught it

Kinetic energy from the same shell integral reproduces the history to **1.0000 on all
three components** and mass to 0.999999, using the code's own stretched faces
(`int r^2 dr = (r_out^3-r_in^3)/3`, `int sin th dth = cos th_l - cos th_r`, and SUM over
longitude -- an early version averaged over longitude and then multiplied by ONE cell's
dphi, under-counting every integral by exactly n_phi = 128). Because KE gated exactly and
ME did not, the discrepancy had to be in the ME definition rather than in the geometry.
**Gate an integral on a quantity the code prints, then trust only what passed.**
Cf. [[validate-the-instrument]], [[sp-hydro-vs-mhd-comparison]].


## AMENDMENT 2026-09-04: bcc cannot even show WHERE B_r is, once the ratio is large

The advice above -- "use bcc only to show WHERE the field is" -- was followed for
sp_dhj_ctl at rot 2 and gave a confident, WRONG answer ("77 % of B_r energy within 5 deg of
the poles"), retracted within the hour ([[sp-polar-field-blowup]]). At rot 2 the
dump/history ratio for the radial component is **0.007** (it was 0.113 at rot 100 above):
bcc holds 0.7 % of the B_r energy, and a quantity holding 0.7 % of a thing cannot locate it.
For B_r on these runs, WHERE requires FACE fields, i.e. a RESTART file, not a dump.

The rule that would have caught it in advance: before localising any energy from a dump,
integrate it over the dump and check it reproduces the history TOTAL. If it does not, the
dump quantity is not the one the history integrates, and the map is of something else.
