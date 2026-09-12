---
name: red-giant-wb-polytropic-walk-overflow
description: ROOT CAUSE of the prod11 3.2106e7 death AND (very likely) the "WB momentum kick at a cold cell" (2026-09-10 05:50, B10 capture): the polytropic well-balanced walk (wb_background.hpp WBAdvance wb_opt==3) extrapolates T linearly in Phi with the stencil's own dT/dPhi; across a 1-cell cold spike (8941 -> 2874 K) tm -> ~0, 1/tm in the exponent, exp overflows to a FINITE 1e191-1e295, WBGuard (non-finite/<=0 only) passes it, the deviation reconstruction does inf-inf -> NaN p and e at the interface; fmax(NaN,pfloor) masks p, HLLC 0*NaN poisons the energy flux
metadata:
  type: project
---
Chain (all captured verbatim by nan_report instrumentation, B3/B7/B8/B9/B10, session 2c3cf987):
cold 1-cell spike (T 2874 K, rho 4x, in a 9000 K layer) -> WB polytropic walk to the ip1 slot:
bg_d 1.2e295, bg_e 8.2e302, bg_p 2.3e302 (finite!) -> deviation w0-bg PLM-limited and re-added
-> dr(IDPR) = NaN, wr(IEN) = NaN at the cell's faces (hydro_fluxes.cpp, before the fmax mask
at :252-253 which turns the NaN p into pfloor=1e-12 and hid this for 2000 cycles / 333 events)
-> hllc_hyd.hpp: flx(IEN) = qc*fl.e + qd*fr.e with qd=0 and fr.e=NaN -> NaN (momentum flux
finite: it never uses w(IEN)) -> u0(IEN) NaN in 3 cells -> rg_grav's eos.Pressure(d, w0(IEN))
NaN -> momentum NaN -> FATAL. Fixes proposed, being implemented 06:00 (agent):
1. wb_background.hpp wb_opt==3: clamp |dlntdphi*dphi| <= 0.5*t0 (smooth fall-back toward the
   isothermal member where the stencil's T gradient is unresolvable; mirrors the isentropic
   sibling's u>0 guard).
2. WBGuard: also reject finite-but-absurd walked states (d,e,p outside 1e-6..1e6 of the anchor).
3. hllc_hyd.hpp general branch: guard wl/wr(IEN) before el/er (mirror the cp clamp).
4. hydro_fluxes.cpp:252-253: explicit !isfinite check instead of the fmax mask (and for IEN).
Gate: V4's first two dumps must stay BIT-IDENTICAL (healthy cells never trip the clamps);
B11 = rg.00124 with the fixed binary, WB ON, no wb_rmax/efloor: death gone? and is the
v_r(305) sawtooth of [[red-giant-prod11-died-3e7]] (the "WB kick") gone too?
Related: [[red-giant-wb-kills-ambient-medium]] (same family), [[red-giant-floor-energy-creation-fix]]
(FL4: wb_rmax alone cured the open-top death = consistent with this being the same bug).

## FIXED AND GATED 2026-09-10 (agent afee698; hunks read and verified by me)
Working tree (uncommitted): wb_background.hpp (a_eff clamp: T may not fall below 0.5*t0 over a
segment; WBGuard absurd() = d,e,p outside [1e-6,1e6] x anchor -> flatten), hllc_hyd.hpp
(eil/eir fallback to p/(G1-1) when w(IEN) not positive-finite), hydro_fluxes.cpp x1 floor
(non-finite p or e -> pfloor and pfloor/(G1-1) consistently). Binary: build_guard/src/athena.
Gate A: build clean, cpplint counts unchanged. Gate B (G1_bitcheck 195376): dumps 00000/00001
BIT-IDENTICAL to V4 (clamps never fire in a healthy star). Gate C (B11_wbfix 195377, rg.00124,
WB ON, nothing else): passes both death cycles, dt 30.14 flat to tlim 3.215e7, ZERO nan_report
lines (B10: 333 wb_recon + 344 hllc), and the v_r(305) SAWTOOTH IS GONE (max consecutive-cycle
jump 3.2 cm/s vs 917 in B3) => the "WB momentum kick at a cold cell" WAS this overflow.
B11_wbfix/bin per-cycle dumps DELETED 05:35 (user ok); the track is in _analysis_0910/task_B11_track_out.txt.
NOTE: this fix does not touch the cold-collapse -> conduction-dt problem (B4/B5); the clamp
opac_tmin=3200 remains necessary for the lidded run.
