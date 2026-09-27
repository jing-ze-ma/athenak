---
name: red-giant-vpert-seed-chart-imprint
description: 2026-09-10 02:00 — the cube-vertex "chart asymmetry" is NOT a seam/vertex code bug: problem/vpert seeds v with a sinusoid in MESHBLOCK-LOCAL (x2,x3), the same pattern stamped into all 96 blocks, so the flow carries the chart's symmetry from t=0 (prod11 3-copy spread 2e-2..0.15 at the FIRST dump; the {0,2},{1,4},{3,5} panel pairing is that seed's exact residual symmetry). With vpert=0 or a chart-free seed all 3 copies of all 8 vertices stay BIT-IDENTICAL with a real flow. Fix = problem/vpert_cart=true (src_vtx, bin_vtx) and re-seed FROM SCRATCH; restarts keep the imprint. V4_prod12cart 195275 tests it
metadata:
  type: project
---

Agent aa3e87412c077b925 (session 293c860c). Reproducers (1 node, 16 ranks, ~660 cycles,
red_giant/V1_sym, V2_seed, V3_cart): worst same-vertex 3-copy rho spread over all 320
radii = 0 (vpert=0), 6.75e-3 (production seed), 0 (chart-free seed with max|v1| 1.2e5 and
tangential flow crossing seams) -> the tangential seam transform, corner fill, seam-end
resample and gnomonic geometry are CLEARED. The asymmetry is a seam-RING phenomenon (seam
edge spread ~130x interior; the vertex sits at the 65-94th percentile of seam cells), i.e.
the vertex is only where two seams meet. In V2 the 4 corners of ONE panel agree to exactly
0.0 until 1.2e6: the state has no stochastic noise, only deterministic symmetry.
Code: src_vtx/src/pgen/red_giant.cpp: vpert_cart (default false, bit-identical), radial
seed v1 = vpert cs sin(3 pi-ish radial) * (4(k4-0.2)+6 k6), k4 = x2y2+y2z2+z2x2, k6 = xyz of
the unit direction (PanelToCart); v2=v3=0. Diff reviewed by me: correct. Build trap: bare
cmake picks gcc 7.5 -> pass -D CMAKE_CXX_COMPILER=$(which mpicxx); test suite needs
module anaconda/3/2023.03. cs test tst/test_suite/rad/test_rad_cs_raddiff_cpu.py passes.
Verification RUNNING: V4_prod12cart (job 195275, 2 nodes, from scratch, full T15/prod12
config, vpert_cart=true, open top, sponge off, tlim 3e6, ~4 h): does the vertex chimney and
the 1.43114e6-type death appear at all? Note the lidded prod11 had the same imprint and no
chimney: the open top is still what makes the seam ring drain; the seed decides WHERE.
Related: [[red-giant-vertex-chimney]], [[red-giant-vertex-floor-hydro-trigger]].

## V4 FINAL 2026-09-10 04:40 (agent report, key lines checked): seed alone is NOT a cure
V4 (old binary: no floor fix, WB active, no wb_rmax, opac_tmin 2500, open, sponge off) died at
9.4456e5, INSIDE T6's 7.7e5-1.04e6 window (T6_nodust's own NaN: 8.62e5): conduction-dt
collapse in a cold-collapsed INTERIOR cell (rank 19, gid 57, panel 3 (J,K)=(10,18), i 277-280,
r 3.39e12), then dfloor + T 3.16e10; dt frozen at 1e-41 from cycle 61400 (job cancelled 04:41).
It bought time: first sub-1000 K cells at 8.8e5 (72 cells, 327 K) vs T6's 4.4e5; clean to 8.4e5.
Vertex copies bit-identical (spread 0.0) at 4e5 and 8e5; 2134/2208 cold cells interior, 2 vertex.
NO floor event ever (e/rho ~1e9 >> 1.42264e7): the open-top death path here is the molecular-
knee cold collapse -> conduction dt, NOT the WB/floor NaN. => prod12 needs the seed AND
wb_rmax + efloor_from_ekin AND the opacity clamp (opac_tmin 3200; from-scratch behaviour of
3200 is the open question). Files: _analysis_0910/task_V4final_*.
