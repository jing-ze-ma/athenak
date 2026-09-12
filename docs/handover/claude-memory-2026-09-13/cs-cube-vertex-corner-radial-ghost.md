---
name: cs-cube-vertex-corner-radial-ghost
description: The cubed-sphere CUBE-VERTEX corner halo (FillPanelCornersFC) only ran over the ACTIVE radial range, so splitting x1 left the radial-ghost layers of that corner block written by nothing -- 0.2 of |B| stale; this was the whole of the "resistivity across a radial block interface" defect
metadata:
  type: project
---

Committed **397b4ad3**, 2026-08-30. Closes the item [[cubed-sphere-resistivity]] left open.

## The defect

`MeshBoundaryValuesFC::FillPanelCornersFC` (`src/bvals/bvals_fc.cpp`) extrapolates the
ng x ng corner ghost block at a CUBE VERTEX, because the generic 3D corner buffer there is
non-reciprocal and its exchange is skipped (`IsCubeVertexCorner`). Its `par_for` looped
**i over `is..ie+1`** -- the active radial range only. The cell-centred twin
`FillPanelCornersCC` has always looped `0..n1-1`, the whole array.

With ONE MeshBlock spanning x1 that is invisible: x1 ends on a physical boundary at both
ends, the corner slot has no neighbour at all, and the physical BC writes those cells
afterwards. **Split x1 and the (radial ghost) x (cube vertex) block is a real ghost region
that NOTHING writes** -- so it kept stale values. One-line fix: loop the full i range.

## Why it hid, and why resistivity found it

**Ideal MHD never reads that block** -- which is exactly why a radial split reproduced
ideal L1(B) to seven digits and the three ideal cubed-sphere tests are BYTE-FOR-BYTE
unchanged by the fix. The **gnomonic resistive curl does** read it: `BcovXi`/`BcovEta` in
`resistivity_gnomonic.cpp` sample at `(i-1, j-1)` and `(i-1, k-1)`.

## Numbers (iprob=11, EVOLVED FIELD gate, changing only `<meshblock>/nx1`)

|  | Linf(B) | heating ratio |
|---|---|---|
| nx2=16, 1 block | 5.8850e-04 | 1.025324 |
| nx2=16, 2 blocks BEFORE | 7.4405e-04 | 1.026503 |
| nx2=16, 2 blocks AFTER | **5.8850e-04** | **1.025324** |
| nx2=32, 1 / 2 / 4 blocks AFTER | 2.3112e-04 all | 1.006358 all |

Bitwise at 1/2/3/6 MPI ranks. The corner-halo Linf itself went **2.09e-01 -> 1.71e-03**
(815 face values were more than 1e-2 off the exact static solution).

## The instrument: cs_test's GHOST SCAN

New in the same commit, and the thing that found this in one run. B is static, so every
HALO cell has a known exact value; the scan reports Linf of that error split by WHICH
directions are ghost -- interior / radial / tangential / **r x EDGE** / **r x CORNER** --
because each names different machinery (x1 face buffer, panel-seam halo, x1x2+x3x1 edge
buffers, cube-vertex corner). MPI-reduced. **This is the only gate in cs_test that sees a
halo region no interior norm reaches.** Use it first for any cubed-sphere boundary bug.

**METHOD, and it worked:** the previous session's localisation (Linf on x2f at
i = 0, j = 0, k = ke) named the RING where the radial interface meets the tangential block
edge. The step that closed it was asking *what fills the ghost cells that stencil reads*,
and then measuring the halo directly instead of the evolved field. Also decisive: an
`eta = 1e-10` run showed the corruption is **resistivity-independent**, which turned it
from an EMF question into a boundary-values question.

## GPU: PASSES (job 11247242, apudev/vipa1001, 2x MI300A, xnack+)

`bench/cs_resist_gpu/`. Every figure matches CPU to all printed digits: 1 vs 2 radial
blocks both Linf 5.8850e-04 / ratio 1.025324 on one GPU, 2 GPUs identical, and the ideal
`mhd_ref_16` (2.931187e-04 / 3.158684e-03) and `mhd_ctl_16` (5.277894e-04) are unchanged
from their pre-fix recorded GPU values. Forced a rebuild of both changed files so warnings
could not be hidden by an incremental build; clean for gfx942.

## THE GATES THEMSELVES WERE BROKEN TOO (b4db75d1)

Found by the GPU/MPI check, not by any run failing.

* **Rank-local**: the EMF L1/max, the panel-edge-ring split, the Ohmic heating ratio,
  max|v| and the div B maxima were all computed over `nmb_thispack` and printed
  UNREDUCED. Two ranks reported heating 1.011584 and 1.039064 for the run whose true
  value is 1.025324. Now reduced before dividing; identical at 1/2/3/6 ranks.
* **L1(B) double-counted the shared x1f layer**: 39936 faces on one radial block vs 41472
  on two, the 1536 difference being exactly the 6 x 16 x 16 duplicated x1 faces. That
  alone moved L1 1.4292e-05 -> 1.3847e-05 for an IDENTICAL solution. Now skipped on the
  block with a same-level inner x1 neighbour, and **L1 agrees exactly across the split**,
  which turned 397b4ad3's arithmetic argument into a measurement.

## Refinement + resistivity: GUARD LIFTED, c4735f18

Measured with the refusal temporarily disabled, iprob=11 with a radial level boundary at
r = 1.5, against the unrefined control:

| nx2 | refined Linf / ratio | control Linf / ratio |
|---|---|---|
| 16 | 5.8397e-04 / 1.019298 | 5.8850e-04 / 1.025324 |
| 32 | 2.3016e-04 / 1.004858 | 2.3112e-04 / 1.006358 |

i.e. it converges at the control's rate, lands just below it, and the heating ratio
approaches 1 faster than unrefined -- the same evidence pattern that lifted the ideal MHD
refinement guard. Telescoping is at round-off (max|dM| 5.3e-23 over 1536 coarse/fine faces, scale 3.2e-07),
the seam flux likewise (dM/dt 2.6e-19), the run is BITWISE over all five dumps at 1/2/3/6
MPI ranks and on 1 and 2 MI300A GPUs, and all four remaining guards still fire. New
`inputs/tests/cubed_sphere_resist_smr.athinput`; iprob=11's final check now also runs
CSTestConsSums / SeamFlux / LevelFlux, which it never had. The level-boundary seam halo is
still FIRST order and is tracked separately in
[[cs-crosslevel-seam-halo-first-order]] -- it is NOT specific to resistivity.

## Still open

`mhd_fofc.cpp` has no cubed-sphere support (latent, `mhd/fofc` defaults false); viscosity,
conduction and AMR likewise, all FATAL at startup. Resistivity + REFINEMENT is still
REFUSED in `Mesh`'s ctor; its stated reason (the radial interface) is now closed, and the
guard's comment still cites the retracted pre-exchange `efld_resist` numbers.
