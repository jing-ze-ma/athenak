---
name: red-giant-prod11-launched
description: "PRODUCTION prod11 launched 2026-09-09 ~06:40 (chain 194515 -> 194516 -> 194517, 6 nodes, 12 h each): the 1.1 R domain on the proven 320-cell grid with every fix of the session -- seed only above 0.5 R (vpert_rmin), ck_nquad 2, no RT clipping, direct RT source, no MLT, no ambient medium. Pinned binary prod11/bin_prod11/athena. The old-physics chain (prod_topre) was cancelled, including its running job 194263."
metadata:
  type: project
---

/orion/ptmp/jinma/Athenak/red_giant/prod11/  (fresh start, tlim 1e8, bin every 2e5, rst
every 5e6; sub.sh restarts from the newest rst so resubmitting it continues the chain).

Config = prod_topre/rg.athinput (1.1 R, x1max 3.52e12, nx1 320, c = -0.395781 +9.763785
-16.467666 +7.899662, measured 5.25/5.27/5.40/5.80 cells/H_p at i=0/2/5/10, median 12.3)
plus: vpert_rmin 1.6e12, ck_nquad 2, rt_de_max -1 (rt_src_direct/rt_newton default on),
mlt_alpha 0, no bg (bg_rho default 0, wb_rmax default 0).  eos_logd_min -14, dfloor
1e-18 unchanged (min rho ~1e-11 on this domain).

Why C (user's choice 2026-09-09): the 1.1 R star is stable and physically sound with the
session's fixes; every remaining failure came from the 1.5 R + cold ambient medium, whose
free fall piles up and crushes the thin atmosphere ([[red-giant-15R-grid-starved-deep]]).
Option A (hot hydrostatic radiatively decoupled corona, T ~ 6e5 K, rho matched ~5e-22,
kappa = 0 above the join, conduction cut off there) is the proposed follow-on -- NOT yet
approved by the user ("hot" changes the spec) and not started.

Risk to watch: at 8x8 the seed-cut/nquad-2 star (R2a) died at t = 1.44e6 in a transonic
thin atmosphere (Mach 2.5 at r/R 1.074); production at 32x32 with the old physics never
showed it.  If prod11 dies near 1.4e6 that pathology is not an 8x8 artifact.
What to check first: out.txt startup lines (nquad 2, vpert_rmin, direct source), dt ~30 s
flat, L_rad,out/L in the face-budget prints (expect 0.3-0.6 early), deep v_rms/v_mlt
staying < 0.05 inside 0.3 R, surface v_rms/v_mlt reaching 0.5 by ~1e6.

**UNCOMMITTED working tree (2026-09-09 07:00)** -- everything the session changed is in
`git diff`, nothing committed (the user has not asked): src/hydro/{hydro.hpp,hydro.cpp,
hydro_fluxes.cpp} (<hydro>/wb_rmax), src/utils/two_stream_rt.hpp (direct source in all
three kernels + blend handover, rt_semi_lin, rt_newton, rt_src_direct, rt_face_flux
accessors, srcraw debug field), src/pgen/red_giant.cpp (bg_rho/bg_temp/bg_rtop ambient
medium, opac_floor + table extension to log rho -24, guard refuses only the dense side,
accretion outer ghost in the background, vpert_rmin, mlt_rmax/tau taper/mlt_mean/
mlt_ramp_time/mlt_x_thr/mlt_relax_time and the chi-limiter removal in the mean path
[should be reverted if mlt_mean is revived], L_rad,out face-budget print, rt_apply_debug
wiring, vc_d_ release), tools/solar_convection/eos_dump.cpp (grid-range args).  prod11
runs on a PINNED copy of the binary built from this tree.  Before committing: revert the
mean-path chi-limiter removal, and consider dropping mlt_rmax/tau-taper (superseded).

## Report 1 (t = 4.4e5, 2026-09-09 07:15): healthy
dt 30.65 flat (drift 2e-4), zero warnings; mass exact to 7 digits; L_rad,out/L settled at
0.31-0.32 (L_cut/L ~0.10); DEEP QUIET PASS: v_rms/v_mlt = 0.0000 at i = 5-20, 0.0003 at
40, 0.017 at 60 -- the seed cut holds at 32x32.  sup at i=5..8 is O(2e-4) but STATIC (dt
~1e-6 per 2e5 s; it was already O(2e-4) at t=0 in this analysis, unlike the old growing
dipole) -- treat as an IC/np.gradient offset at the RCB kink; watch.  Surface: the seed is
decaying (i=180: 0.31 -> 0.07), no convection yet at 4e5 (grey_prod had it rising by
6e5-9e5).  Photosphere 3780 -> 3407 -> 3458 K at 1.068-1.072 R (dipped 10 %, recovering).
Report 2 due at t >= 1.5e6 (crossing the 8x8 transonic-death time 1.44e6).

## R6_open (job 194526, 07:25): open outer boundary test at 8x8
User's suggestion.  R2a config (1.1 R, seed cut, nquad 2, no MLT) + problem/outer_bc = open
(Stein-Nordlund continuation, outflow allowed, inflow clamped: open_outer_noinflow = true),
tlim 2e6.  Question: R2a died against the WALL at 1.44e6 (transonic thin atmosphere,
Mach 2.5 at r/R 1.074) -- does letting gas leave fix that, at what mass-loss rate, and does
the star stay sane?  Pinned binary at red_giant/bin_R6_open/athena (the first submission,
194524, failed instantly because that path did not exist yet).  Monitoring agent reports
at 1.6e6 or death together with prod11's Report 2.  prod11 confirmed nx2 = nx3 = 32.

## Option A APPROVED as a test (user, 07:40): R7_hot_open, being implemented by an agent
Hot hydrostatic corona: above the join (star's column at 2e-22) continue the hydrostatic
integration ISOTHERMALLY at bg_temp = 6e5 K from p_join (pressure-matched; density drops
by T_star/T_c -> light on top, RT-stable; H = kT r^2/(mu m_H GM) ~ 2 r so it sits still).
New: problem/bg_hydrostatic (default false), problem/opac_corona (0) + a conduction/RT
cutoff radius <hydro>/rad_kappa_rmax (default 0) so the corona neither radiates (t_cool at
kappa 1e-5 would be 0.4 s) nor sets the conduction dt (K ~ 1/kappa).  With the direct RT
source, kappa = 0 gives e0 = 0 -> src = 0 exactly; the Newton branch is skipped (Em = 0).
Outer boundary: problem/outer_bc = open (user's call) -- Stein-Nordlund continuation,
inflow clamped.  Test on the 480-cell 1.5 R grid at 8x8, tlim 1.5e6.  Pass = corona
|<v_r>| < 1e5 and T within 20 % of 6e5, dt >= 15 s, photosphere within 10 % of 4000 K,
deep quiet.  wb_rmax stays at the join; dfloor 1e-23; Bondi at 6e5 K: rho_max ~2.6e-17
vs ~5e-22 pressure-matched.

## R6_open RESULT (08:00): the open outer boundary REMOVES the 1.44e6 death
Ran to tlim 1.5e6 with zero collapses (R2a, the identical walled star, died at 1.447e6).
Both are transonic (Mach 2.1-2.3 in the 20 cells above tau=2/3, 8e5-1.4e6) -- the open BC
does not remove the supersonic 8x8 surface flow, it survives it (the flow leaves instead of
reflecting and shocking).  Photosphere (3757 K at 1.067 R at 1.5e6; peak 4932 K at 9e5),
L_out/L (0.34 -> 0.66 -> 0.38) and deep quiet (0.002-0.004) identical to R2a; net mass lost
4.6e27 g = 5.2e-5 Msun/yr averaged (a burst tied to the surface transient; 1e-7 of the
envelope; shuts off by 1.5e6).  => outer_bc = open is the safer production default; if
prod11 (wall) shows the dt dip deepening near 1.4e6, restart its next chain segment with
outer_bc = open (in the input file; restarts embed the input, so it needs a fresh input --
check how -r handles a changed parameter before relying on it).

## Report 2 (t = 1.58e6, 08:40): PASSED the 1.44e6 point
dt 30.6507 (drift 2.6e-4 from cycle 0), zero warnings; mass unchanged to 7 digits.  The
atmosphere above tau=2/3 is SUBSONIC in v_r at 32x32 (Mach 0.8-1.2; max|v_r| 7.3e5 vs
1.2-1.3e6 at 8x8) -> the wall is never struck transonically: the 8x8 death was resolution.
Deep quiet saturated (i<=20: 1e-4; i=60 plateau 0.029 after 1e6); RCB values frozen to 4
digits (static IC offset).  Surface convection developing first: v/v_mlt at i=210/240/270
0.10/0.09/0.07 (4e5) -> 0.58/0.65/0.47 (1.4e6), i=150 decaying.  L_rad,out/L 0.31 -> 0.66
by 1.56e6.  Photosphere 3780 -> 3407 -> 4832 (1e6) -> 5207 K (1.4e6) at r/R 1.077 -- hot
vs T_eff 4000: read with the LID in mind ([[red-giant-grey-opaque-lid-bug]]).  Decision on
a lid-free fresh restart waits for R8_lid_open/R8_lid_wall (jobs 194540/194542).

## 09:50 decision: prod11 is superseded by R9 if R9 passes
prod11 = lidded + walled (stable only because of the lid).  Its chain continuations
194516/194517 CANCELLED; 194515 runs on as the lidded control / fallback until R9's
1.5e6 report.  If R9 passes: cancel 194515 and start production FRESH from the R9 config
(tlim 1e8 + restart chain) -- R9 cannot be extended (tlim 2e6 < first rst at 5e6).  If R9
fails: resubmit prod11's chain (sbatch prod11/sub.sh with --dependency=afterany:194515).

## 10:00: R9_lid_open32 DIED at t = 1.99e5 (hydro dt 37 -> 7e-15 in one step, rank 23 =
panel 1, then FATAL non-finite).  A NEW 32x32 failure: the same lid-free open star lived
to 1.3e6 at 8x8.  prod11 (lidded, walled) is again the only stable production run: its
chain was RE-QUEUED (two 12 h continuations after 194515).  prod12 (R9 config, tlim 1e8)
is staged but NOT to be submitted until R9's death is understood.  Monitoring agent is
localizing the cell (open-ghost at the top? photosphere? seam/vertex? deep?) from the
t = 1e5 dump and comparing prod11 at the same cell/time.
R9 death signature (from its nan_check prints): 1,636,536 non-finite cells at cycle 6500
with "bad i span ~200..323" on EVERY rank, 26 cycles after a one-step hydro collapse
(37 s -> 7e-15) with dt ~1e-15 in between -- nothing can move; only the two-stream's
DOWNWARD sweep from the TOP ghost spreads a NaN to every i >= icut (~200 here) in one call
on all panels.  Suspect: the OPEN ghost fill_open_out's DensFromPT(p_g, T_i) failing
(1e6 sentinel/NaN) for a 32x32 top-cell state (8x8 open survived to 1.3e6).  Fixes being
implemented: (1) validated open ghost with fallback to the column ghost + a fire counter;
(2) the RT top ghost takes rho/T/kappa/B from the top ACTIVE cell ie, not the hydro ghost.
L_rad,out/L in R9 before death: 0.58 (3e4) -> 0.36-0.42 (lid-free; prod11 lidded 0.43).
