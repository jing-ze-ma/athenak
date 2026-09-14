---
name: cs-mhd-blast
description: b283ad3a -- the iprob=12 blast gained a uniform tilted field built from a VECTOR POTENTIAL, so div B is zero to round-off by construction; it is the test that actually exercises the fc seam halo under a shock
metadata:
  type: project
---

**b283ad3a.** `cs_test` iprob=12 (the blast) now has an MHD form: a uniform tilted field
built from the **vector potential** `A = 0.5*(B0 x r)`. Each face value is the
circulation of A around **that face's own edges** divided by its area, using the very
same `dxedge` and `area` arrays `mhd_ct.cpp` consumes.

**Why a vector potential and not a face projection.** CT preserves whatever divergence it
is handed, so a field projected face-by-face (what iprob=8 does) starts with a small but
**PERMANENT monopole** on a curvilinear grid. Built as a discrete curl instead, div B is
zero to **ROUND-OFF by construction**: measured 1.5e-15 at t=0, and CT holds it at
<= 1.6e-13 through the shock. Follows `pole_mhd_blast.cpp`.

`A.rhat = 0.5*(B0 x r).rhat == 0` identically, so the RADIAL edge potential A1 vanishes
and only the two tangential edge components are needed -- the same simplification the
spherical-polar version uses.

**Gate the construction with TWO numbers, not one.** div B alone is satisfied by a ZERO
field. Also check the magnitude: `|B|max = 0.09997` against `b0c = 0.1` confirms the
discrete curl reproduces the INTENDED field and not merely some divergence-free one.

**Why it matters:** unlike the radial split monopole the generic branch supplies, this
field has finite **TANGENTIAL** components on every panel, so it is the first test that
exercises the face-centred seam halo, the corner EMF and the CT update **when a shock
arrives**. Survives the panel interior, a seam and the cube VERTEX with 0 non-finite
cells of 49152, on CPU, MPI (np=4, div B identical to serial) and GPU.

Parameters: `problem/b0c` (default 0.1) and `bhx/bhy/bhz` (0.37/0.61/0.70, a deliberately
non-symmetric tilt). Also adds a profile/scatter diagnostic behind the `CS_SCATTER`
env var and the div B gate itself. See [[cs-blast-vs-cartesian]],
[[cs-shocks-through-seams]].
