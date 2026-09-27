---
name: cs-seam-ghost-eint
description: RETRACTED as a mechanism -- the seam degrades EVERY quantity by the same ~15-19x over a block face; eint is not special and the 1/beta signature was in the block-face row too. The E-vs-B halo inconsistency is real but ~14 %. Instrument <mhd>/cs_seam_diag (81a50486), fix b9efbbfa (single-rank only)
metadata:
  type: project
---

**Built and measured 2026-09-04, commit 81a50486.** `<mhd>/cs_seam_diag = N` prints, every
N cycles, how far a SEAM ghost cell's state departs from a quadratic extrapolation of the
three active cells inside it. Lives in `src/mhd/mhd_seam_diag.cpp`, called at the end of
`MHD::ConToPrim` where w0 and bcc0 have just been filled over the ghost zones.

## The mechanism, and why it is not just "the resample is inaccurate"

The cell-centred halo carries the CONSERVED variables (`SendU` packs `u0`, whose IEN is the
TOTAL energy with `0.5|B_src|^2` inside it) while the face field crosses the seam through a
SEPARATE exchange with its own co-location, gnomonic transform, along-seam resample and
monotone clamp. ConsToPrim then forms `e_int = E_ghost - KE - 0.5|B_ghost|^2` from two halves
resampled by DIFFERENT operators -- the energy as a scalar, the field componentwise and then
squared. Resample-then-square != square-then-resample, and the whole residual lands in e_int,
divided by beta.

## The numbers (cubed_sphere_mhd_strat, iprob 13, 2x2 blocks per panel, cycle 1)

    region        ncell      dens mean        eint mean        magE mean
    SEAM  beta<1  19200      7.85e-07         2.23e-04         1.17e-05
    SEAM  beta>1   5376      2.11e-06         7.52e-06         6.95e-06
    block beta<1  19200      4.66e-08         1.18e-05         6.94e-07
    block beta>1   5376      1.47e-07         4.45e-07         2.58e-07

**Read the block-face row correctly**: at a same-panel block face the ghost is a plain COPY,
so that row is NOT halo error -- it is the INSTRUMENT'S OWN extrapolation error, i.e. the
null control. Subtracting it, the seam low-beta eint anomaly is ~2.1e-4 against ~7.4e-7 for
density in the SAME cells: **a factor of 285 between two quantities that went through the
same resample.**

**The beta split is what makes it a mechanism rather than a correlation.** Across beta = 1
the seam eint anomaly grows 30x (7.5e-6 -> 2.2e-4) while density's moves the OTHER WAY
(2.1e-6 -> 7.9e-7). A resampling error hits both alike and does not care about beta; a
residual from subtracting a separately-resampled 0.5|B|^2 is amplified by exactly 1/beta.

## Why this is now the leading candidate

It is the only surviving mechanism that fits ALL of: boundary-associated and seam-hottest
([[cs-ulp-amplification]]), worst at low beta, immune to every interior accuracy fix (the
well-balanced source, GS07 EMF, PPM all leave the HALO untouched), cured only by dissipation
(which damps the resulting error rather than removing it), and a floor trigger IN THE GHOST
CELLS -- which is why `eos_efloor` overflows before death and why the counters cover ghost
zones. It also survives the rotation fix, which removed a different confound entirely
([[cs-rotation-source-bug]]).

Found by a Fable 5.1 agent reading the halo code; the numbers above are this session's own
measurement of its hypothesis, not the agent's claim.

## THE REMEDY IS BUILT -- and it is only 14 % of the anomaly (b9efbbfa)

`<mhd>/cs_seam_econsist` exchanges the field-free part of the energy, `E - 0.5|B_cc|^2`,
across the seam and rebuilds each seam ghost cell's total energy from it plus the ghost's
OWN field. Measured on the same test:

    seam beta<1  eint   2.23e-04 -> 1.92e-04   (14 % removed)
    seam beta>1  eint   7.52e-06 -> 5.39e-06   (28 % removed)
    dens, magE, and BOTH block-face rows: unchanged to three digits

**So the mechanism is REAL and SUBDOMINANT.** The null control passes -- the correction
moves only what it should -- but ~86 % of the seam low-beta anomaly is something else, still
unidentified. **Downgrade the "leading candidate" claim below: it is a confirmed real defect
that does not by itself explain the instability.**

A hypothesis tested and REFUTED: that the residual was the along-seam MONOTONE CLAMP acting
differently on E and on 0.5|B|^2, two large numbers whose difference is small -- if so,
resampling the small quantity directly would fix it. It does not: both formulations give
1.92e-04 to three digits. The committed version resamples the small quantity anyway, being
the better conditioned of two equal-cost forms.

**NOT MPI-SAFE.** One rank is validated end to end (flag off is BIT-IDENTICAL to the
pre-feature binary over 20 cycles). TWO ranks HANG DURING SETUP -- before any stage runs, so
it is in the new MeshBoundaryValues object's construction or its first InitRecv, not in the
correction kernels. UNRESOLVED. A guard refuses the flag above one rank rather than leaving
it switchable, since every production run here is multi-rank and a hang would be blamed on
the physics. **Anyone picking this up: start with that hang.**

Left alone deliberately: the KINETIC term has the same structure (momentum transformed
exactly then resampled, so KE in the ghost is square-of-resample against resample-of-square)
but grows with the Mach number rather than with 1/beta.

## The original remedy sketch, for reference

Send `e_int` (or `E - 0.5|B_cc|^2`) as the seam CC quantity and rebuild
`E_ghost = e_int + KE(ghost) + 0.5|B_cc,ghost|^2` AFTER the FC halo is filled -- a dual-energy
treatment restricted to the halo. One extra pass in `RecvAndUnpackFC`/`FillPanelCorners*`.
Gate it by re-running this diagnostic: the seam low-beta eint row must fall to the density
row's level.

## Traps this instrument has

* It needs **MORE THAN ONE MeshBlock PER PANEL**, or every tangential boundary is a seam and
  the null-control row comes back EMPTY. With nx2=nx3=32 pass `meshblock/nx2=16
  meshblock/nx3=16`. The test input documents this where the flag is defined.
* The block-face row is extrapolation error, not halo error. Do not read it as "the block
  face is also broken".
* beta here is the proxy `2*eint/B^2`, not `2p/B^2` -- it is a bin boundary, not a
  measurement.
* Gated: 20 cycles with the flag on and off are BIT-IDENTICAL, and asking for it off the
  cubed sphere is refused (verified to fire).
