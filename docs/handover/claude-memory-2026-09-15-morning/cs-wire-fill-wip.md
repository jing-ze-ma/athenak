---
name: cs-wire-fill-wip
description: The cube-vertex corner fill: BOTH open problems fixed and now PROMOTED to mesh/cs_vertex_fill, default ON (branch cs-wire-wip, 437e9b43). 3-4x on the corner halo, ideal MHD bitwise unchanged, MPI-bitwise. MERGED to polar-average-perf; widening icoar/ifine MEASURED not worth building. START HERE
metadata:
  type: project
---

**Branch `cs-wire-wip`, HEAD `437e9b43`** (dafd6997 fixed the two problems, 437e9b43
promoted it). The env var is GONE: the switch is **`<mesh>/cs_vertex_fill`, default TRUE**
on a cubed sphere, a no-op on every other grid, and setting it false reproduces the
previous answers exactly. **NOT merged to `polar-average-perf`** (b283ad3a), which is
clean and unaffected. Read the three commit messages; they are long and complete.

**IDEAL MHD AND HYDRO ARE BITWISE UNCHANGED** -- their stencils never read the cube-vertex
corner block. Verified on cubed_sphere_mhd: the dumps differ in exactly 5 bytes, the
`true`/`false` string of the input echoed in the HEADER (a dump carries its input, so
never conclude from `cmp` alone). Only the RESISTIVE curl reads that block, which is why
the resistive inputs are the only recorded numbers that moved.

**A SHOCK IS UNDISTURBED:** the MHD blast at a cube vertex, contrast 100, is identical on
every gate with the fill on or off (0/49152 non-finite, div B round-off both ways);
contrast 1000 survives, and so does the same blast with eta = 1e-3.

## THE DESIGN, unchanged and sound

A cube vertex is where only THREE panels meet, so the ng x ng corner ghost block has no
diagonal neighbour and its exchange is skipped. Widen the buffers that already cross the
seam so the corner rides the existing exchange: no new slot, no new tag, MPI-safe by
construction. `angles()` works in the SOURCE chart and the along-seam resample extends
continuously past the active range, so **no new 2-D machinery is needed**.

## BOTH open problems are now FIXED

1. **OWNERSHIP is EXACT, not FindPanel.** The seam between the two flanking panels is
   *exactly the DIAGONAL* of this panel's own gnomonic chart: for P0 = +z, tan xi = x/z
   and tan eta = y/z, and the +x/+y boundary is the plane x = y, i.e. xi = eta. So with
   `a = s2*xi - pi/4`, `b = s3*eta - pi/4`, a > b belongs to the x2-face panel and a < b
   to the x3-face one. **Complementary by construction** -- which is the whole point: the
   old FindPanel test had a GAP, and because FillPanelCornersFC stops extrapolating once
   the wire is on, a cell both panels rejected was written by NOTHING. That, not
   degeneracy, is why "ownership OFF" measured better. Take each component at its own
   STAGGERING (LeftEdgeX for the axis it lives on) and restrict to a TRUE cube vertex.
2. **The RADIAL GHOST corner** is reached by widening the x1x2 (16-23) and x3x1 (32-39)
   EDGE buffers by the same rule -- a triply-ghost cell can only be reached by a buffer
   already ghost in x1. Their seam normal matches the flanking face's, so the resample
   carries them unchanged. Gate a radial-ghost layer on BOTH edge neighbours being at
   this level; one alone leaves half the block unwritten (the 397b4ad3 trap).

## MEASURED, on the mesh that actually exercises it

`cubed_sphere_resist.athinput`, nx2 = 16, PLAIN RADIAL SPLIT (`meshblock/nx1=4`):

```
              tangent ghost   r x CORNER SEAM (int / phys)   evolved L1 / Linf
  OFF           1.5968e-03      1.1176e-03 / 6.4703e-04      4.0656e-05 / 8.0656e-04
  wire          4.2842e-04      3.7000e-04 / 6.4742e-05      4.0216e-05 / 7.7989e-04
```

On `cubed_sphere_resist_smr` (every radial interface is ALSO a level boundary, so only the
tangential half applies): tangent ghost 1.6455e-03 -> 4.1699e-04, evolved Linf 5.6994e-04
-> 5.4221e-04. The WIP's own best was 1.1193e-03. **Only `isame` is widened, so at a level
boundary the wire correctly does not apply** -- that is the remaining gap, and it needs
`icoar`/`ifine` widened too.

## A NEW GATE, because the old one could not see this

"r x CORNER SEAM" pooled cells whose radial neighbour EXISTS with cells OUTSIDE the domain.
The second set is larger, set the maximum, and MASKED the first -- exactly what the wire
moves. cs_test now reports the two separately and prints the gid/level of the two edge
neighbours at the worst interior cell. See [[validate-the-instrument]].

## VERIFIED

Bitwise dumps at 1, 2, 4 MPI ranks with the wire ON, no hang; the no-op control bitwise
against OFF under MPI; HIP compile clean; style clean.

## DECIDED: promoted, default ON (the user chose this)

**COST:** the x2/x3 face and x1x2/x3x1 edge messages grow by 2*ng along the seam, ~25% on
a 16-cell block and ~6% on a 64-cell one. **GAIN:** 3-4x on the cube-vertex corner halo,
3-5% on the evolved resistive Linf, nothing at all for ideal MHD or hydro.

**ONE PLACE IT IS WORSE:** `cubed_sphere_resist_smr_narrow` (MeshBlocks four cells wide)
loses 1% of evolved Linf, because that mesh falls back to a plain index copy along the
seam so the data the fill carries is itself first order. Recorded in that input.

**WHAT IS LEFT:** widen `icoar`/`ifine` so the fill also applies at a coarse/fine boundary,
and merge the branch to `polar-average-perf`.


## MEASURED 2026-09-02 on HEAD 5c369988: WIDENING icoar/ifine IS NOT WORTH BUILDING

The "WHAT IS LEFT" item above was measured before building it, on
`cubed_sphere_resist_smr` (the mesh where every radial interface is also a level
boundary, so the wire's radial-ghost path is switched off by its own level test).
Bench `/viper/u2/jinma/ATHENAK/bench/cs_buf`, fill ON vs OFF:

    category            nx2=16 OFF    nx2=16 ON   gain      nx2=32 OFF   nx2=32 ON   gain
    tangent ghost        1.7207e-03   2.9512e-04  5.83x     1.3558e-04  7.8982e-05  1.72x
    r x CORNER SEAM      2.1689e-03   2.1698e-03  1.00x     6.1636e-04  6.1652e-04  1.00x
    EVOLVED L1           1.1068e-05   1.1070e-05  1.00x     6.6727e-06  6.6726e-06  1.00x
    EVOLVED Linf         3.5288e-04   3.5288e-04  1.00x     1.2352e-04  1.2352e-04  1.00x

**The argument is the pair of rows, not either one alone.**  `r x CORNER SEAM` is
indeed the largest halo category and indeed gets NOTHING from the wire (1.00x) --
that is the gap, and it is real.  But `tangent ghost` is its sibling, fed to the same
resistive curl, and improving it **5.83x moved the evolved field by 0.00%**.  A
comparable improvement to the corner would therefore buy nothing measurable either.
Halo quality in these categories is simply not the binding term on a refined mesh.

And the corner category is not a stalled defect: 2.1689e-03 -> 6.1636e-04 is a ratio
of 3.52, i.e. converging at ~1.8 order like everything else.

**The 5% the wire used to buy on this input is GONE because the SEAM CO-LOCATION FIX
subsumed it.**  The recorded "evolved Linf 5.6994e-04 -> 5.4221e-04" predates 999458a3;
the same run now reads 3.5288e-04 either way.  So the co-location fix both improved the
number 1.6x and removed the wire's margin on it.

**DECISION: do not widen `icoar`/`ifine`.**  The cost is intricate index code in the
send/recv path where ndat must match across a coarse/fine pair with f1/f2 sub-block
offsets -- the most MPI-critical code in the repo -- for a measured payoff of zero.
Revisit only if a problem appears whose answer depends on the cube-vertex corner halo
at a level boundary.  Cf. [[cs-cube-vertex-real-exchange]], which was declined on the
same grounds, and [[measure-impact-before-claiming]].
