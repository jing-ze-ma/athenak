---
name: lowmach-solver-audit
description: Audit of the all-Mach hydro solvers lhllc (Minoshima+Miyoshi 2021) and ausmpup (Liou 2006 + Edelmann 2021) -- algebra matches the papers; two design choices (LHLLC chi from the NORMAL velocity only; AUSM+-up pressure diffusion hard-wired to Mref=1 so f_a^p = 1, the criticised p-u asymmetry) -- and the rk2-vs-rk3 stability test on the stratified column at rest (2026-09-07): NO CFL restriction found, rk2 == rk3
metadata:
  type: project
---

**Algebra (read against the papers, 09-07 06:00):** LHLLC's contact pressure expands to
p* = [m_R p_L - m_L p_R + phi m_L m_R (u_R - u_L)]/(m_R - m_L), the M&M form, phi=1 gives
Batten's HLLC; S_M untouched (correct). AUSM+-up: a*, a-tilde, M-bar, f_a, alpha, beta,
M4+-, P5+-, K_p/K_u/sigma, p_u x P5+P5-, upwinded mass flux -- all Liou 2006. Design
choices, not bugs: (1) LHLLC chi = max|u_n|/max c uses the NORMAL velocity (commented
draft used the magnitude) -> in shear the fix is over-applied; verify vs the paper.
(2) AUSM+-up M_p uses f_a^p with Mref = 1 hard-wired (Edelmann 2021) => f_a^p == 1, the
1/f_a of Liou is gone: pressure diffusion at high-Mach strength, velocity diffusion ~M --
exactly the p-u asymmetry critics cite. Liou's balance could be restored now that the WB
deviation reconstruction makes face dp a perturbation; Mref should be an input.

**"CFL must scale with M" objection, TESTED** (wb_column + problem/vpert seed, commit
after ebd57244; scratch lowmach/; stably stratified tgrad=-0.1, WB polytropic on,
etotgrav on, Mach 1e-3 seed with a grid-scale (-1)^(i+j+k) part, 400 time units ~ 235
box crossings, hst KE + odd-even projection of dumps):
- 3D CFL 0.3: hllc damps KE 1e10x; lhllc and ausmpup 4e3x then hold; rk2 == rk3 to 3
  digits for both; the grid-scale mode persists ~100 units under the fixed solvers
  (94 % of KE at t=100, 11 % at 400) but never grows.
- 2D (nx3=1) CFL 0.45 (legal 2D limit 0.5): same, rk2 == rk3. Mach 1e-4 seed: same.
- 3D CFL 0.5 and 0.8 kill hllc too (unsplit 3D limit is 1/3) -> void, not a solver thing.
- TRAP: the stock wb_column input has tgrad = -0.3, SUPERADIABATIC (ad = -0.286): every
  solver "grows" KE 300x by convection. Use tgrad = -0.1 for a rest test.
So: on this grid no CFL restriction beyond the ordinary one, with rk2. My von Neumann
estimate (CFL <~ M^(1/3) for rk2) was pessimistic: the PLM and the remaining dp/(rho c)
coupling supply the missing damping. Not yet tested: Gresho vortex, sound-wave damping
rate, and these solvers in a dhj run.

**UNSEEDED column (user's suggestion, vpert = 0, same setup, to 3000 units):** the flow
that develops is purely VERTICAL (v2 = v3 = 0 exactly, so no checkerboard possible) and
concentrated at the bottom wall: it is the WB scheme's discrete residual driving a
steady column flow, and the low-dissipation solvers let it grow to a larger steady
amplitude instead of damping it. KE saturates: hllc 1.2e-16 (v_rms 5e-8, Mach ~1e-8),
lhllc 1.5e-13 (Mach 2e-6), ausmpup 5.7e-12 (Mach 1e-5, rising 8e-13 -> 5.7e-12 over
t = 100 -> 1600 then FLAT). rk2 == rk3 throughout. No instability; the cost of the fix
is a 1e2-1e3x larger residual velocity floor at the wall (AUSM+-up 6x LHLLC).

## The two papers (read 09-07 06:30-07:30), and what they mean for us

**Velasco-Romero & Teyssier 2025 (MNRAS 537, 2387; arXiv 2405.11063) Appendix B:**
2D isothermal atmosphere (rho = p = exp(-x), g = 1, gamma 5/3, [0,1]^2, N x N, periodic y,
fixed hydrostatic ghosts in x), FV2 = MUSCL-Hancock + their LHLLC + their WB scheme. HLLC +
WB holds the profile; LHLLC + WB GROWS an instability, growth rate rising with C and with
resolution; they derive from a modified-equation model with the WHOLE upwind dissipation
scaled by M and a first-order-in-time anti-diffusion -Delta t/2 d_tt u: nu ~ |a| h (M - C)
=> stable only for C < M ("impossibly stringent"). Their Gresho at M = 1e-3 needed C = 0.05.
OUR TEST of the mechanism: same atmosphere, N = 64, our LHLLC, WB isothermal, rk1 (forward
Euler) vs rk2 at C = 0.3: rk1 blows up for HLLC too (KE 1e-25 -> 1e-5 by t = 40; LHLLC
faster, 1e-6 by t = 10); rk2: both flat at round-off for 60 crossings, N = 64/128/256,
rk2 == rk3. So their C < M is the anti-diffusion of a (effectively) first-order-in-time
integrator, which SSP-RK2 does not have; and their model over-counts the fix (LHLLC
scales only the Delta-u pressure term, the Delta-p mass-flux term stays O(1)). Not a
property of LHLLC + RK2.

**Edelmann+2021 (A&A 652, A53) Sec. 4 rest tests:** isothermal and isentropic
atmospheres, rho-p linear reconstruction, 64 cells, C = 0.9, mostly implicit ESDIRK23,
run to 5000 t_BV (~16000 t_SC). 1D: every WB method holds Mach < 1e-12, flux choice
irrelevant. 2D: all WB runs start at 1e-13 but the LOW-MACH AUSM+-up ones grow
EXPONENTIALLY over a few thousand t_BV (pattern: resolved horizontally, grid-level
vertically; "pressure-velocity decoupling + gravity => an instability very similar to
convection, but in stable stratifications"; also with Miczek+2015 and Li & Gu fluxes);
the standard AUSM+B-up does not. Their remedy is pressure diffusion (Edwards & Liou),
i.e. the p-term, with an INDEPENDENT cut-off Mcut^p ~ 0.1 so 1/f_a^p stays bounded
(their eq. 38) -- our Mref = 1 is the strongest version of that (f_a^p = 1). Their
timescale: e-fold ~150 t_BV; for our isothermal test t_BV = 10, t_SC = 0.77 => needs
1e4+ time units: isoL_N64_{hllc,lhllc,ausmpup}_rk2 + lhllc_rk3 launched 07:30 to
t = 2e4 (2000 t_BV), ~9 h each on one core, scratch lowmach/isoL_*.

## AUSM+-up cut-offs are now INPUTS (hydro/ausm_mcut, ausm_mcut_p; defaults 1e-13, 1.0)
Carried in EOS_Data (the solver signature is fixed). Tested 09-07 08:00, seeded
stratified column (Mach 1e-3), unseeded column, isothermal 2D atmosphere:
- Liou ORIGINAL (shared cut-off 1e-3, M_p ~ 1/f_a ~ 500x): NaN within ~10 steps on
  EVERY test, even unseeded. Shared 0.1 (5.3x): NaN with any seed (smooth or grid-scale,
  Mach 1e-3 or 1e-4) in ~7 steps, fine unseeded; Edelmann 0.1 (mcut 1e-13, mcut_p 0.1):
  the same. mcut_p 0.3 and 0.5: fine to t = 400. mcut_p = 1 (default): fine.
- MECHANISM (diag_ed1, dumps every step): the bottom WALL cell drains -- rho 0.91 -> 0.12,
  eint to the floor in ONE step at t = 0.032-0.039. At a mirrored wall M4+(M)+M4-(-M) = 0
  so the wall mass flux is rho a M_p alone: AUSM+-up's pressure-diffusion term makes a
  solid wall LEAK mass in proportion to the cell's pressure perturbation with
  coefficient K_p/f_a^p; at 5x it is fatal within steps. Edelmann's "stability issues
  when 1/f_a^p becomes exceedingly large" is this. So: the original AUSM+-up is
  UNUSABLE with walls in a stratified atmosphere here; the decoupled cut-off is
  required and 0.1 is on the edge; keep mcut_p >= 0.3 (default 1).
- Unseeded floors: liou1 (0.1 shared) is the QUIETEST (KE 1e-15 column, 1e-27
  isothermal, = hllc) because the stronger p-term damps the residual; ed1 ~ default.

**RESOLUTION (09-07 09:00, commit "AUSM+-up: cut-off Mach numbers and the wall-HLLC
hybrid as inputs"):** the mcut_p = 0.1 death was NOT the wall (it recurs with
hydro/ausm_wall_hllc = true, HLLC on the physical wall faces, which is correctly keyed
on mb_bcs = the mesh's flag only at the domain edge). It is a TIME-STEP limit: M_p is an
explicit pressure diffusion with coefficient K_p/f_a^p, stability ~ CFL < f_a^p/(2K_p):
0.1 -> 0.38 (marginal at 0.3: dies; FINE at CFL 0.1 and 0.05); 0.3 -> 1; 1 -> 2. Liou
original with M_inf = 1e-3 dies even at CFL 0.002. So "CFL ~ M" is TRUE for the
original AUSM+-up (via the p-term) and Edelmann's decoupled cut-off is precisely the
cure; FALSE for LHLLC + RK2 (no such term). The bottom cell fails first but the wall is
where it shows, not why. Keep ausm_mcut_p >= 0.3 (default 1).

## WAVE TESTS (09-07 09:30, scratch lw/ and lowmach/gw_*)
- Linear sound wave (built-in linear_wave, 3D oblique 64x32x32, 5 periods, rk2, C 0.3,
  build_generic): total L1 error hllc 1.96e-4, lhllc 1.64e-4, ausmpup 1.57e-4 at amp 1e-3;
  IDENTICAL ratios at amp 1e-6 (errors scale linearly) -> the fixes make sound waves
  slightly LESS dissipated, no phase pathology, no amplitude dependence down to Mach 1e-6.
- Gravity waves (isothermal 2D N=64 atmosphere, smooth Mach 1e-3 seed, 200 units = 20
  t_BV): after the transient, hllc keeps Mach 2e-5 that is 99 % gravity-wave-band power
  and 98 % smooth (k<=4 vertically). lhllc/ausmpup(1)/ausmpup(0.3) keep Mach 1.8e-4,
  STEADY from t=50 to 200, of which ~85 % is NEAR-ZERO frequency (not waves), drho/rho
  1e-4 with only 35 % of its vertical variance in k<=4 and 13-27 % of vx in the vertical
  ODD-EVEN mode: Edelmann's "resolved horizontally, grid-level vertically" decoupled
  pattern, present with BOTH low-Mach solvers, AUSM's p-term strength (1 vs 0.3) makes
  no difference. Not growing over 20 t_BV; the 2000-t_BV isoL runs answer growth.
