---
name: red-giant-sponge-cs-kinetic-energy-bug
description: 2026-09-09 — REAL bug, fixed 887c241e; sponge-off SURVIVES (T12) but the vertex excursion is started by the GHOST-FILL basis bug (T14 testing); T11 with the fix died the same way at 8.72e5 — the red_giant sponge computed KE with the orthogonal sum on the cubed sphere's covariant momenta, draining internal energy at the panel corners; the 24 cube-vertex cells inside the sponge hit the T floor at t~1.0e6 once the open top went transonic -> holes -> NaN (T6/T7/T8/T9)
metadata:
  type: project
---

Chain: T8 (nan_check every cycle) caught the birth: 72 active + 64 ghost cells at i=310-311,
all in the panel-corner blocks (lx2,lx3 in {0,3}), first active cell at the block corner
(k,j)=(2,2)/(2,9)/(9,2)/(9,9) = the cell touching a CUBE VERTEX. T9 fine dumps: at 1.021e6
the vertex cell at i=302 is at e/rho=1.42e7 (T floor) with neighbours at 1.8e12; the hole
climbs one cell per ~3e3 s (304 @1.027e6 ... 308 @1.039e6) while v_h ramps 3 -> 39 km/s;
rho 40-300x below neighbours; neighbours fall in at 3-5 km/s. No anomaly at all before
1.0e6 (vertex v_h <= bulk, zero cold cells). i=302 is the first cell inside the sponge
(zs = rin + 0.98 (rout-rin) = 3.4515e12, i~298).

Cause (src/pgen/red_giant.cpp rg_sponge): ei = E - 0.5 d (v1^2+v2^2+v3^2) with v from u0
momenta -- but on cs u0(IM2,IM3) are COVARIANT components on a non-orthogonal basis; the
right KE is 0.5[m1^2 + (m2^2+m3^2-2c m2 m3)/(1-c^2)]/d with c = pcoord->cos_cell (history.cpp
151-158). The wrong KE drains the true e_int by (KE_true-KE_orth)(1-fac^2) per step: largest
at panel corners, ~v^2, lethal once the top is transonic (8 km/s at i>=300). prod11's vertex
cells are 5-8% colder than the bulk in the sponge layer -- same defect, sub-lethal there.

Fix: CsKinetic() helper + sponge uses it, scales momenta by fac, E = ei + fac^2 ke (exact
on any basis). Built 887c241e 08:45; test T11_spongefix (restart from T9's 1.0e6 rst,
must pass 1.04e6 with no cold vertex cells) handed to the Opus agent. Other orthogonal-KE
sites in the pgen (IC 1460, inner-boundary passes 2229-2314, open ghost 2666) use tiny
velocities or w0 velocities of unknown basis -- audit before production, not urgent.
The agent's sponge positivity guard (skip damping when ei<=0) is a compatible symptom fix.

**How to apply:** on the cubed sphere never form KE from u0(IM2,IM3) with the plain sum; use
cos_cell. Any cs pgen with a velocity sponge/damping that edits u0(IEN) has this bug.
Related: [[red-giant-nan-check-hides-origin]], [[cs-nonorthogonal-audit]],
[[cs-mhd-c2p-floor-corrupts-ue]].

## UPDATE 18:40 — the fix is right but does NOT cure the vertex collapse
T11_spongefix (887c241e, restart from 1.0e6... actually from T9's 1.0e6 rst; the agent's A/B
used matched dumps) reproduces T7 to 3-4 digits in the vertex columns (the fix changes the
vertex |v_h| growth by NOTHING measurable) and died at 8.72e5 (chaotic timing) with the
identical two-dump event: vertex column (m0,k7,j0) i=310-314: |v_h| 3e4 (8.2e5) -> 7.3e5
(8.4e5, e/rho still 1.17e12) -> 3.15e6 with e/rho at the floor (8.6e5), all 8 cube-vertex
columns simultaneously and to IDENTICAL values (the setup is symmetric under the cube group,
so the 8 vertices are exact copies -- one deterministic mechanism). The bulk column in the same
block is smooth. So the sponge's KE error was a real but sub-dominant drain; the vertex
instability lives elsewhere: something that acts on the cube-vertex ACTIVE cell in the
transonic layer. Candidates: (1) FillPanelCornersCC's extrapolated corner ghosts being read
after all (transverse radiative-conduction flux? cs seam flux correction? the along-seam
resample's clamped window at the seam END?) -- DECISIVE CHEAP TEST: poison the corner ghosts
(1e30) in a 1-cycle job and see whether any ACTIVE cell changes; (2) the seam-end resample
quality (bvals_fc comment: halo error grows 7.5x toward the vertex); (3) a sponge-region
effect: rerun with problem/sponge=false from the 8e5 rst (T11 has rst every 1e5) -- if the
vertex still collapses, the sponge is fully exonerated; (4) the two BUG-LIVE ghost fills
(only i=322/323, probably not it, but fix anyway). The 'T = 3.16e10 K' at i=293 is a
DOWNSTREAM solver saturation (10^10.5 bracket), not the origin.

## UPDATE 2026-09-09 ~20:00 — discriminating runs from T11's 8e5 rst (all in red_giant/)
- T12_spongeoff (195050, problem/sponge=false): SURVIVES to 9e5 (T11 died 8.72e5). The vertex
  |v_h| excursion at i=300-320 still happens (peak 5.2e5 at 8.4e5, 2x T11's) but RELAXES
  (8.6e4 by 9e5). Sponge = the killer of an excursion it does not start.
- T13_vfilloff (mesh/cs_vertex_fill=false): NULL, dies at the identical 8.721147e5 --
  FillPanelCornersCC is called unconditionally for CC (bvals_cc.cpp:688-691); the flag only
  widens FC buffers.
- Poison test (build_poison, <mesh>/cs_corner_poison=true, 2 cycles vs control, bit-
  reproducible): 232 active cells change, all 24 vertex corners, max rel 2e-6, via eint
  (conduction cross term conduction.cpp:384-393 reads (ks-1,js-1), limited). NOT the driver.
- Code study: the only O(1) vertex-only defect left is the pgen's radial ghost fills
  (fill_open_out/fill_open/fill/IC) writing contravariant v into the covariant u0(IM2,IM3)
  and an orthonormal KE; at c=-0.451 ConToPrim returns v_h x1.82 and e_int - 0.82 rho v_h^2
  in the ghost (= audit's BUG-LIVE 1+2). Fix = GnomonicEquiangleLowerMom formulas with
  c = cs ? cos_cell(m,k,j) : 0. Test T14_ghostfix (tree copy /orion/ptmp/jinma/Athenak/src_vfix,
  build_vfix) running from the 8e5 rst, tlim 1.2e6.
