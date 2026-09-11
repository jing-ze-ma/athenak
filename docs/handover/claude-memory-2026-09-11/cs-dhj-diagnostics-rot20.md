---
name: cs-dhj-diagnostics-rot20
description: cs MHD dhj VALIDATED at rot 20 against cs hydro and sp hydro -- 3D T structure identical on all grids, equatorial jet weak (~1 km/s at 1e-3..1e-2 bar) on EVERY grid, the field does NOT slow it (beta >= 100 everywhere); report URL + the analysis pipeline (EOS inversion, cs panel->lat/lon, contravariant wind projection)
metadata:
  type: project
---

2026-09-06 evening. Report (artifact): https://claude.ai/code/artifact/1a0f4a2b-f35d-42e0-b376-ed938ac3e76a
Runs compared at rot 20: cs_prod_hyd_rot (dump 00010), cs_inflow_bcfix (00020, fixed bottom
BC 248f1b77), sp_dhj_hyd (00010, old binary); time series to rot 60 from cs_prod_hyd_rot,
cs_prod_mhd_rot (drifting BC, rot 0-26) and sp_dhj_hyd.

**Results.** (1) T on the 1e-3/1e-2/0.1/1 bar isobars and T(p) at substellar/antistellar/
terminators are IDENTICAL cs-hydro vs cs-MHD and match sp: no panel/vertex imprint. A
nightside equatorial hot spot at 0.1 bar, lon ~+140 deg, is on all three grids (physical,
unexplained). (2) Zonal-mean zonal wind: narrow eastward equatorial jet 0.5-1.1 km/s
between 1e-4 and 1e-2 bar (oscillates rot to rot), broad westward mid-lat flow above
1e-3 bar, interior < 0.3 km/s below 0.1 bar. The day-night flow is 15 km/s at the top but
divergent. THIS SETUP HAS ONLY A WEAK JET on every grid (3.5 d rotation, UHJ). (3) B does
NOT slow the jet through rot 20: MHD jet inside the hydro scatter, identical below 0.1
bar; beta >= 100 at 0.01 bar, > 1e4 deeper; |B| 0.1-30 G at 0.01 bar (weakest under the
substellar point), 10-300 G deep. (4) max|u| at 0.01-0.1 bar is 25 % lower on cs than sp
= the 3x coarser equatorial cell (2.8 vs 0.84 deg); the zonal mean is unaffected.

**Pipeline ($S/dhjcs.py, figs.py, fig3.py, build_page.py; $S = scratchpad of session
967431f9).** EOS: `<mhd>/eos_table_dump` on a serial CPU run of the sp input with tlim=0
(build_dhj_cpu, mesh 32x8x32/mb 32x8x16, ~1 min) gives eos_table.txt; the inverter masks the
NaN rows (T < 71 K) and was validated to 5e-4 in T and p against problem/ck_dump_file
(needs time/nlim=1, ck_dump_j/k set). cs geometry: MB m is on panel m//(nmb/6) (gid order
follows panel_trees); xi = pi/4 * x2, eta = pi/4 * x3 with x2,x3 in [-1,1]; unit vector
(a x + b y + n)/delta with cubed_sphere.hpp PanelFrame; substellar is -x, lon =
atan2(-y,-x). w0 vely/velz are components on the UNIT tangent vectors e_xi, e_eta
(non-orthogonal, cos = cos_cell up to 0.475): V = v1 rhat + v2 e_xi + v3 e_eta. bcc is
(B_r, B.e_xi, s B^eta) on the orthonormal pair (e_xi, e_perp). sp lat needs
StretchTheta: theta = pi/2 (1 + sinh(a(2 xi-1))/sinh a), xi = x2/pi, a = 3.
TRAPS: isobar interpolation must CLAMP its weight to [0,1] (non-monotonic p in floored
top columns extrapolated to 7000 km/s -- no cell exceeds 21 km/s); 5-deg lat/lon bins are
empty near the poles (nearest fill). See [[cs-mhd-bottom-inflow]], [[sp-hydro-vs-mhd-comparison]].
