---
name: cs-test-ffdecay-rotaxis
description: cs_test iprob 14 (force-free resistive decay, twin of sp_test iprob 11) and rot_axis variants for cs rigid rotation (4990eb41); first cs-vs-sp tables and the reporting traps
metadata:
  type: project
---

Committed 4990eb41 on `polar-average-perf` (pushed to fork), 2026-09-06.

**What exists.** `src/pgen/cs_test.cpp` iprob 14 = the exact decaying force-free field
B = b0 (sin az, cos az, 0), faces from GL4 Stokes loops of A = B/alpha along every panel
edge (div B ~ 1e-14), user radial ghosts at the decayed state, finalizer `CSTestFFCheck`
(L1/Linf v, p, face B by panel interior / seam / cube vertex, div B, resistive EMF on
every edge vs eta*alpha*B; appends `<basename>-errs.dat`). Input
`inputs/tests/cubed_sphere_ffdecay.athinput`. Rigid rotation (iprob 3) takes
`problem/rot_axis = z | edge | vertex | custom` (vertex = axis through a cube vertex, the
pole's analogue); `inputs/tests/cubed_sphere_rigidrot.athinput` now sets `rot_axis` and
`conv_errors = 1` (without conv_errors the by-region L1 table is NOT printed at all).

**Results (uniform, defaults, n = cells per panel edge, nx1 = n/2, t = 1 rot / t = 2 ff):**
- rigid rotation global v rates: z 2.2/2.4, edge 2.0/2.1, vertex 2.0/2.1; the cube-vertex
  region 1.6-2.2. Versus sp: rigid z 2.6, rigid x 2.3 with polar rows 1.5. cs has NO
  1.5-order region; the vertex is no worse than the sp pole.
- force-free decay (eta 0.05): global v 1.8, p 2.0, B 2.2 at n8->16; seam == interior;
  EMF error 2-3 % at n8, ~2nd order. sp: v 1.9, B 2.4, polar v 1.1. n32 pending at
  hand-off (`$S/csff/n32`, ~70 min on a loaded node).

**Toroidal (iprob 11, bunif 0, `cubed_sphere_toroidal.athinput`) and blast, added later
the same day (fefe6a17 +):** toroidal ideal B 1.94/1.95 (seam 1.9, vertex 1.6->1.2),
resistive with user ghosts 2.1 and Ohmic ratio 1.017/1.006/1.002 (sp 0.90/0.98/0.99);
mass round-off on reflecting walls. TRAP: b0c must keep p0 - b0c^2 R^2 > 0; b0c 0.5 in
r <= 2 puts p = 0 on the outer equator and the outer layer reads first order. Blast, same
problem as sp (open shell, `cubed_sphere_blast_open.athinput`): panel centre vs edge /
vertex 0.5 / 0.4 % L1 at n16, 0.08 / 0.15 % at n32; cs centre vs sp equator 0.4-1 %, cs
vertex vs sp pole 0.15-0.5 %. The cube vertex costs about what the sp pole costs. The
cs-vs-sp blast comparison needs the same tlim (hydro input default is 0.30, sp 0.2).

**Traps.** A command-line override (`problem/rot_axis=...`) needs the parameter to EXIST
in the input file or parameter_input.cpp fatals at line 390. iprob 3 on cs also prints
CONSERVED SUMS with Lx/Ly/Lz (equal for the vertex axis, a free consistency check).

**Why:** the user wants the same simple exact-solution matrix on cs as on sp so the two
grids can be compared on equal footing before any dhj claim. See [[sp-test-rigidrot-resist]],
[[sp-pole-fixes]], [[cs-hydro-validation]].
**How to apply:** run `$S/csmat.sh`-style series with nx1 = n/2 and all three doubled
together; read tables with the collector pattern in cscollect.py (interior/seam/vertex).
