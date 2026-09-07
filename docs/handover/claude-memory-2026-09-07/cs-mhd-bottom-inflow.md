---
name: cs-mhd-bottom-inflow
description: GATE PASSED 09-06 (bcfix 0 inflow cells through rot 20, ctl 2 at rot 19-20 as predicted; cs_mhd_prod launched from scratch). cs_prod_mhd_rot died at rot 32 by a bottom-boundary inflow; CAUSE FOUND AND FIXED (248f1b77): the dhj radial BC built the cs ghost bcc as plain face averages while ConsToPrim subtracts the orthonormal-frame energy, so the ghost energy DRIFTED every call. Every cs MHD dhj run before 248f1b77 carries it
metadata:
  type: project
---

2026-09-06. cs_prod_mhd_rot went NaN at rot 32.1 (cycle ~642000, 5.4 h into job 11470115)
and ran on with NaNs to rot 57 (rst 00004-00006 = rot 36/46/56 are NaN; 00003 = rot 26
is finite but inside the eruption; 00001 = rot 6, 00002 = rot 16). Jobs 11470115/6
cancelled by the user's instruction at 07:35.

**What it was.** NOT the cube vertex: 86-89 % of the radial KE is in panel interiors, 0 %
in the vertex band. A radial INFLOW at the inner boundary cell (v_r = -2.1 km/s at i = 0,
decaying to -0.9 by i = 6) in 2 isolated columns at rot 19 (0 at rot 9, -0.1 at rot 15),
at lat +-10-15 deg on the equatorial panels 0/1/2/4 at all longitudes; by rot 29 a band.
The columns carry a 400 G flux tube anchored at the bottom (B_r -388, B_xi -275 G at
i = 0-1) and lose 20 % density per 10 rot; global mass -1 % by rot 30 (sp -0.2 %), KE1 20x
sp, |B| max 5 kG at i = 0 before the NaN. The neighbours of an inflow column have v_r = 0.

**sp control.** sp_mhd_diss at rot 20 has |B_r| 100-136 G at i = 0 in the same band (cs
23-30 G) and NO inflow anywhere (min v_r -0.02 km/s at i = 0 through rot 50). The bottom
field grows to kG on BOTH grids: that is the setup, not the defect.

**Differences left.** equatorial resolution (sp theta-stretched 0.84 deg vs cs 2.8 deg),
the cs non-orthogonal MHD flux/EMF path at the boundary face, `<mhd>/cs_lowbeta_llf`
(default ON on cs only), the cached cs well-balanced source. The pgen BC
(HydrostaticEquilibrium, deep_hot_jupiter_rt.cpp 2406-2972) IS complete: inner ghosts =
the initial hydrostatic column with zero momentum (restarts carry ghost zones), ghost faces
mirrored about the boundary face (tangential with a sign flip, flux-scaled), ghost energy
corrected for the field each call. An earlier claim that it wrote no ghost field was WRONG.
The only defect found was the IC zeroing w0(IVX) three times (fixed after 18f5dd21;
harmless). The boundary is therefore the same on both grids.

**Submitted 08:05: bench/cs_inflow_{ctl,hlle,nollf}, jobs 11502664/5/6**, 979edada binary,
restart from rot 6 to tlim 20 rot, bins every rot, 6 h apu1: ctl must reproduce the 2
inflow cells by rot 19 (gate = `count of bottom cells with v_r < -1 km/s`, script pattern
in this session); hlle and nollf say whether the Riemann solver / fallback matters.
See [[cs-mhd-dhj-blowup]] (the vertex mechanism list), [[sp-mhd-energy-excess]].

**Experiment results so far (09-06 ~10:30).** hlle: inflow ALREADY at rot 11 (36 cells,
97 by rot 14), at HIGH latitude (panel 3 = north polar panel), after the deep field
(peak i = 16-24, 400 G) migrated to the bottom cells between rot 9 and 10 (<ME> at i = 2-4
2700 -> 11400 erg/cm3 by rot 11, 65 G mean B_tan at i = 0). ctl and nollf are BITWISE
identical through rot 12 (the low-beta fallback never fires at the bottom). Bottom cells
are at REST on both grids (|v_h| 0.01 km/s at i = 0-5), so no wall shear layer: the field
reaches the bottom from above and piles up because the wall is Neumann in B_r (sp too:
|B_r| 114 G at i = 0, no inflow). More dissipation -> earlier/larger pile-up and inflow.
Submitted: bench/sp_theta_uniform (job 11504808, sp_mhd_prod input with f_stretch_theta
0.001 = uniform 2.8 deg theta, the cs cell; from scratch, 20 rot) and bench/cs_inflow_cart
(job 11504809, cs_cart_momentum = true, restart rot 6).

**CAUSE FOUND (09-06 09:00), FIXED 248f1b77.** The user asked whether a covariant/contravariant
transform was missing at the cs boundary -- it was, for the FIELD. The pgen BC rebuilt the
ghost bcc as plain face averages (components on the non-orthogonal tangent pair) while
ConsToPrim/GnomonicEquiangleRaiseVelMHD stores bcc in the orthonormal frame and subtracts
0.5|bcc|^2 there; the inner-ghost energy update does `E -= 0.5|bcc_old|^2` (orthonormal,
left by C2P) then `E += 0.5|bcc_new|^2` (plain) -> the ghost TOTAL energy drifts by the
cross term (c^2(b2^2+b3^2) + 2c b2 b3)/s^2 EVERY CALL. Verified from the restart files
(they store ghost zones; parser pattern: header = par_end+10+232+20, lloc 24 B/MB, per-MB
[u0(5,20,20,132), x1f, x2f, x3f]): ghost rho spread 0.0 but ghost e_int spread 1e-3 (rot 6),
8 % (rot 16), 3.4e8..2.2e9 about 1.46e9 with negatives (rot 26). Low-pressure ghosts drain
the bottom cell. sp has c = 0; hydro has no field. Fix: rotate the ghost bcc into the
orthonormal frame in both ghost kernels. HIP rebuilt (make -j 6, 09:22). Experiment
bench/cs_inflow_bcfix (restart rot 6, new binary) submitted; ctl/hlle/nollf keep running
as negative controls; cs_inflow_cart and sp_theta_uniform CANCELLED (superseded).
CONSEQUENCE: every cs MHD dhj result before 248f1b77 (cs_prod_mhd_rot, rcmfix, the
sheet/eruption archive, the 0.2-rot vertex deaths at nx 64?) ran with a drifting bottom
boundary; the production cs MHD run must be redone FROM SCRATCH on >= 248f1b77.
The INITIAL CONDITION and the restart bcc rebuild had the same plain-average form: fixed in
the commit after 248f1b77 (IC energy off by 1.9e-3 of e_int in the top cells at t = 0,
zero at depth; one-off, not a drift). A from-scratch cs MHD run needs the HIP binary built
after that commit.
**CPU VERIFICATION DONE (10:40):** 300 cycles (4197 s) from the rot-16 restart, old vs fixed
boundary, ghost zones from the written rst: old binary changes the ghost (E-KE) by up to
8.4e5 erg/cm3 at i=0 (rms 3.8e4) and 1.1e6 at i=1 while the real magnetic-energy change is
2.4e3 -- the drift, 300x the physics; fixed binary: max 2.9e3 = the magnetic-energy change
plus round-off. The mechanism is CONFIRMED by direct measurement, independent of the GPU
experiments.

**GPU GATE (09-06 19:50, $S/inflow_gate.py = count of i=0 cells with v_r < -1 km/s per bin):**
ctl 0 through rot 17, 1 at rot 18, 2 at rot 19-20 (min v_r -2.4 km/s, rho_min(i=0) 6.14e-4
-> 4.81e-4, mass -2.3e-3, E -1.5 %); nollf == ctl bitwise (fallback irrelevant); hlle 36
cells at rot 11 -> 367 at rot 20; **bcfix: 0 cells at every rot 8-20, min v_r -0.006 km/s,
rho_min(i=0) flat 6.11e-4, mass -9e-4, E -0.9 %** (the remaining drift is the wall's
O(B^2 h^2) leak + Ohmic, no drain). Bottom |B|max grows the same on both (~300 G at rot
20): the pile-up is the setup, the drain was the ghost energy. CLOSED.
**PRODUCTION: bench/cs_mhd_prod launched 09-06 19:5x** on build_dhj_gpu 09:32 (md5
efd92190.., 551cd78a; cccab476 only touches fofc/refinement): lead job 11516277 (15:20 h,
fitted before the 09-07 12:00 -> 09-12 12:00 maintenance), then 11516227 -> 11516228
(24 h chain). Input = cs_prod_mhd_rot's; bins 2/rot... every 2 rot, rst every 10 rot.
