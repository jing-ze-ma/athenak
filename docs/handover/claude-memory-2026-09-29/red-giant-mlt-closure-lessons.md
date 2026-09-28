---
name: red-giant-mlt-closure-lessons
description: "Per-column per-face MLT (problem/mlt_alpha) is noise-dominated deep (x~1e-5 carries L, gradient noise +-5e-6 -> F flickers 0/L on adjacent faces) and no tau-taper can hand over gently at the top (tau falls 1.1 dex per cell: 0.17 L per cell at best). The fix in progress: SHELL-AVERAGED 1D MLT with f = min(F_MLT(x), L - F_rad,diff)."
metadata:
  type: project
---

2026-09-09, runs R1_mlt / R1b_mlt / R1c_mlt in /orion/ptmp/jinma/Athenak/red_giant/ (8x8
per panel, 1.1 R domain, mlt_alpha = 3 = mlt_alpha_ic), each analysed by an Opus agent.

| run | top treatment | result |
|---|---|---|
| R1 | none (cap 10 % e/step only) | F_mlt runs into the thin layers; photosphere +43-140 %, dt collapse 6.8e5 |
| R1b | hard `mlt_rmax` = 3.42e12 (tau~30) + cap L | 0.79 -> 0.03 L across TWO cells where F_rad = 6e-4 L; T_eff 3763 -> 5808 K; NaN at the cut, t = 6.2e5 |
| R1c | log-tau taper 30..1e4 | drop spread over 12 cells but still 0.17 L/cell at i=283: tau falls 1.1 dex PER CELL there, a 0.05 L/cell taper would need a 22-dex window |

- **Deep flicker.** At i~8-30 (tau~1e12) the x = grad - grad_ad that carries L through
  F ~ x^1.5 is ~1e-5; the discrete two-cell gradient carries +-5e-6 of noise, so F_used
  alternates 1.0 / 0.0 F_req on adjacent faces (drop/L = +-1).  A per-column local
  closure cannot work at this x.  This is also why MLT-on left the RCB pile-up dipole
  exactly as in R0 to 6e5 s (sup(6): 6.1e-5, 1.2e-4, 1.5e-4 at 2,4,6e5 -- marginally
  decelerating, so a working MLT WOULD saturate it; a flickering one does not).
- **Flux budget with MLT on** (R1b, t = 6e5): closed in the radiative core (F_rad 0.97 L)
  and in the CZ i=60-180 (F_mlt+F_conv+F_rad = 0.9-1.1, F_rad < 1e-6 L); collapses to
  0.11-0.19 at i=200-290 = the pile-up under the shutoff.  R0 (no MLT) carries 0.03-0.15
  at i=240-280 -- the resolution-limited surface deficit ([[red-giant-flux-deficit-is-spinup]]).
- **Right design** (being implemented as `problem/mlt_mean = true`): evaluate MLT on the
  shell-averaged profile (noise / sqrt(6144)), apply the 1D F_mlt(r) to every column, and
  close with f = min(F_MLT(x), L/(4 pi r^2) - F_rad,diff(mean)) so the handover happens
  wherever the MEAN gradient lets radiation carry the flux -- no tau window, no radius.
  Per-cell safety cap (10 % of e per step) kept.  This is the 1D-stellar-structure MLT
  plus resolved 3D convection on top; where the resolved flow flattens the mean gradient,
  x and F_mlt drop -- the handover in the mean.

Why the IC is consistent only with THIS: the column was built with the same F_MLT(x) and
F_rad(grad) so F_rad + F_mlt = L at every point at t = 0 (radiative for tau < ic_tau_rad).
Related: [[red-giant-deep-onset-is-rcb-pileup]].

## Shell-averaged closure (mlt_mean) -- R1d, 2026-09-09 03:00

Implemented and stable to t = 2e6 (dt back to 58.9 s after a dip to 21 at 1e6).  The
"photosphere +100-146 % hotter than R0 at fixed i" reading was a TAU-SCALE SHIFT, not a
hot photosphere: the outer layers puffed up, icut moved 281 -> 296, and tau(i=290) went
from 0.055 to 5.5e3, so fixed-i comparisons crossed 4 decades of optical depth.  On the
tau scale (rtcol_late.txt, one column, t = 2e6):

| | T at tau=2/3 | emergent F_lw/F_req | RE at tau=280 would be |
|---|---|---|---|
| R0 (no MLT) | 3045 K (i=287) | **0.198** | -- |
| R1d (mean MLT) | **3794 K** (i=299) | **0.492**, still rising outward | 9600 K; R1d has 8465 (slightly under) |

R0 at 8x8 radiates 0.20 L (0.33 at 32x32): the deficient one.  R1d's cut region sits at
the radiative-equilibrium temperature for carrying L (9911 K at tau~50 vs RE 9930) and
its thin layers are net emitters (cooling toward the 3364 K skin) -- a transient in the
right direction.  Deep faces still drop 0.3-0.7 F_req where x_mean jitters (radial IC
discretization noise, not horizontal); the RCB dipole is unchanged because the closure
trusts F_rad_diff(mean) = 0.99 L at i=6 while the conduction operator delivers ~0.97 L.
Next: L_rad,out over ALL columns as a diagnostic, F_rad from conduction's own discrete
flux, and an mlt_ramp_time switch-on (agent running R1f/R0f).
