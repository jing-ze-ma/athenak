# proto2: 1-D implicit two-moment (E, F) grey radiation prototype, SDIRK2 + ray closure

Files: `scheme.py` (the scheme), `tests.py` (drivers), `*.out` (raw results quoted below).
Python 3.9, numpy/scipy.sparse (`spsolve`), numba for the ray sweeps. CPU, viper login node.

## Scheme as implemented

Equations (c = 1): dE/dt + dF/dx = c ka (B - E); dF/dt + c^2 d(DE)/dx = -c chi F, chi = ka + ks.
Optional gas: rho_cv dT/dt = -c ka (B(T) - E), B = B0 + Bp (T - T0) (linearised a T^4).

Face f between cells L, R (h = stage step = gamma dt, theta_f = 1/(1 + h c chi_f), chi_f mean of the
two cells, a_f = c sqrt(max(D_L, D_R))):

    F_E,f = (F_L + F_R)/2 - alpha_f (a_f/2) [E]_f + w_rc RC_f
    P_f   = c^2 (D_L E_L + D_R E_R)/2 - alpha_f (a_f/2) [F]_f
    RC_f  = theta_f h c^2 [ (g_L + g_R)/2 - ((DE)_R - (DE)_L)/dx ],  g_i = ((DE)_{i+1} - (DE)_{i-1})/(2dx)
    [u]_f = u_R - u_L (1st order)  or  u_R^- - u_L^+ (limited linear reconstruction, minmod or MC)

w_rc = 1 - alpha_f by default (variants: 1 - theta_f, 1). Sources implicit (E: -c ka E + c ka B;
F: -c chi F). All terms are assembled as one sparse 2n x 2n (3n with gas) matrix per pass
(vectorised COO assembly), `spsolve`.

Thick limit (RC check). In a stage, F_i^new = theta_i (F~_i - h c^2 g_i) with F~ the stage's explicit
part. Momentum interpolation to the face, F_f = theta_f (F~_f - h c^2 ((DE)_R - (DE)_L)/dx), equals
(F_L + F_R)/2 + RC_f (with F~_f = mean of the cells). For h c chi >> 1, theta_f h c^2 -> c/chi_f and
F~ theta -> 0, so with alpha -> 0 the E flux -> -(c/chi_f) ((DE)_R - (DE)_L)/dx: the compact
3-point diffusion flux. Numerical check (`thick.out`, first lines): alpha = 0, D = 1/3, Gaussian at
chi dx = 100 vs a separate compact 3-point diffusion code with the same SDIRK2 steps:
max|E - E_compact| = 1.4e-6 (N 64), 2.9e-7 (N 128) (max E 0.7), i.e. the difference is the O(1/tau^2)
transport correction, not a stencil difference.

Time: SDIRK2 (Alexander, gamma = 1 - 1/sqrt 2, L-stable, stiffly accurate). Stage 2 uses
f(Y1) = (Y1 - U^n)/(gamma dt) (exact, conservative). Each stage = npic (default 3) Picard passes; each
pass recomputes the rays (D, C1) and the limiter selection from the current iterate (lagged), then one
linear solve. With reconstruction the limiter is linearised by freezing its selection (which
difference, or 0) from the iterate: the slope is implicit, only the choice is lagged.

Rays: steady formal solution on the current state, 16 Gauss-Legendre mu (8 per hemisphere),
S = (ka B + ks E)/chi (S = 0 where chi = 0), exact integration across cells with piecewise-constant S,
cell-averaged intensities for J, H, K. Periodic: exact wrap (pass 1 gives I_end and the transmission A,
I_0 = I_end/(1 - A)); inflow: beam in the largest mu (I = Eg/w), vacuum on the right.
D = K/J, C1 = sum w |I - S|/J, s_ray = smoothstep((C1 - 0.1)/0.4) (C1-smoothstep, max of the two cells
at a face). Cells the rays leave empty (J < 1e-3 E) take neighbour values (8 passes); beyond that
D = max(1/3, f^2) (1 for a beam, 1/3 isotropic), C1 = 1.

Boundaries: periodic, or left inflow (characteristic ghost: incoming u+ = Eg + Fg/ag given by a beam
consistent with the rays, Eg = 1, Fg = c mu_max, Dg = mu_max^2; outgoing u- = E0 - F0/ag taken from
cell 0) and right outflow (copy).

Options tested (names in the tables):
- alpha: `A1` = (1 - e^-tau)/tau; `A2` = 1 - (1 - A1)(1 - theta); `A2max` = max(A1, theta);
  `A3` = max(A1, s_ray); `A3b` = max(A2, s_ray). tau = chi_f dx.
- `1st` / `minmod` / `mc`: reconstruction of the states in the dissipation terms.
- `rc=oma|omt|one`: RC weight 1 - alpha, 1 - theta, 1.
- `/Dg`: realizability floor D >= f^2, f = |F|/(cE) of the iterate (see the positivity section: NOT
  rays-only).
- `/fb`: step fallback: if a stage is not realizable, redo the step with backward Euler, 1st order,
  same alpha; if still not realizable, BE with alpha = 1, a = c.
- `/a=c`: dissipation speed a_f = c instead of c sqrt(max D).
- `S1` (coordinator addition, linear): alpha = max(A1 psi, s_ray), psi = Jameson sensor
  |u+ - 2u + u-|/(|u+ - u| + |u - u-| + eps) of the stage-start state on E and F/c, face = max of the
  two cells; 1st-order states; ONE solve per stage (`np1`: rays, D, psi from the stage start).
- `rayflux` (coordinator addition, linear): alpha = clip((Fc - Fray)/((a/2)(E_R - E_L)), 0, A1),
  Fc = (F_L + F_R)/2 of the stage-start state, Fray = c sum w mu I_face (the rays' kinetic face flux,
  from the stage-start rays); alpha = 0 where |E_R - E_L| < 1e-12 (E_L + E_R) or the ratio is negative;
  faces whose two cells are both ray-empty take A1. No s_ray. One solve per stage. Runs with the
  piecewise-LINEAR ray source (`--slin`: face source = mean of the two cells, exact integration of a
  source linear in tau); with piecewise-constant S the kinetic face flux of a thick cell is
  -(c/4)(E_R - E_L) (upwind half-range flux, O(tau) too large), alpha clips to A1 and the thick limit
  is the A1 one (kappa_num/kappa - 1 = 0.866 at N 64 and 128, measured).

## Results (all measured; raw lines in the named .out files; finished 16:30 PST 11-14)

Tests: (1) realizability, N 100: beam into a uniform medium, square pulse (F = cE, periodic), two
counter-streaming pulses, beam into a slab (vacuum then chi dx); chi dx 0, 0.01, 0.1, 1, 10, 100,
scattering and absorbing (B = 0); c dt/dx 0.5, 4, 100 (+1000 in `pos_new.out`, `pos_rf.out`); every
stage and step checked, tolerance 1e-12 of the reference energy 1; NaN counts as a violation.
(2) thin wave E = 1 + 0.3 sin, F = 0.2 c sin(+0.5), ka = ks = chi/2, B = 1, t = 0.5, c dt/dx 0.5,
N 64/128/256 vs N 2048 (block averages), error = L1(E) + L1(F). (3) Gaussian (s0 = 0.05) in pure
scattering chi = 6400 (chi dx 100/50/25), t = s0^2/(2 kappa), exact diffusion solution, 16 (N/64)
steps (c dt/dx 96) or 16 (N/64)^2; kappa_num/kappa - 1 from the second moment. (4) checkerboard
0.1 on E = 1, chi dx 100, N 64, 20 steps. (5) gas-coupled linear wave, k = 2 pi, ka = 1 (propagating
mode) and 100 (slowest diffusive mode), rho_cv = 1, B = 1 + 4 (T - 1), D fixed 1/3 (the analytic
dispersion needs constant D) and with rays. (6) sums of E, F (vacuum) and E + rho_cv T, periodic.

### Main table

| variant | test 1 viol. (cases) | worst minE / max(abs F - cE) | thin wave orders chi dx 0.01 / 0.1 / 1 | thick L1 orders, kappa_num/kappa - 1 (N 64,128,256) | time order (succ. diff.) thin / thick | odd-even (cfl 0.5 / 100, 20 steps) | conservation |
|---|---|---|---|---|---|---|---|
| A1/1st | 87 (132) | -0.61 / 0.85 | 1.00 1.08 / 0.97 1.06 / 1.13 1.11 | -0.01 0.00; +0.866 +0.866 +0.866 | 2.02 2.01 2.00 / - | 0.079 / 7e-17 | <1e-14 |
| A1/minmod | 89 (132) | -0.62 / 0.88 | 1.95 1.97 / 1.93 1.98 / 1.99 2.02 | 1.09 1.09; +0.131 +0.065 +0.032 | 2.01 2.01 2.00 / - | 0.079 / 3e-17 | <1e-14 |
| A1/mc | - | - | 1.98 2.01 / 1.97 2.01 / 2.01 2.02 | 2.21 2.20; +0.003 +0.000 +0.000 | - / 1.99 1.79 0.84 | 0.079 / 6e-17 | - |
| A3/1st/Dg | 55 (132); 73 (176 incl. cfl 1000) | -0.31 / 0.69 | - | same as A1 (s_ray = 0 in thick) | - | - | - |
| A3/minmod/Dg | 67 (132); 87 (176) | -0.72 / 1.39 | 1.95 1.97 / 1.93 1.98 / 1.99 2.02 | 1.09 1.09; +0.131 +0.065 +0.032 | - | 0.079 / 3e-17 | - |
| A3/mc/Dg | 94 (176) | -1.03 / 1.94 | 1.98 2.01 / 1.97 2.01 / 2.01 2.02 | 2.21 2.20; +0.003 +0.000 +0.000 | 2.01 2.01 2.00 / - | 0.079 / 6e-17 | - |
| A3/vl/Dg (van Leer) | 105 (176) | -0.87 / 1.69 | 2.01 2.01 / 1.95 2.00 / 2.01 2.02 | 1.95 2.04; +0.029 +0.007 +0.002 | 2.01 2.01 2.00 / 1.99 1.79 0.85 | 0.079 / 3e-17 | <5e-15 |
| A2/minmod | 84 (132) | -0.63 / 0.89 | - | 1.09 1.09 (dt~dx), 1.03 1.00 (dt~dx^2); +0.136 +0.068 +0.034 | - | 0.038 / 8e-17 | - |
| any /Dg/fb (BE fallback) | 42-54 (132); 51-60 (176) | -0.036 / 0.146 | (fallback never fires in tests 2-3) | | | | |
| S1 (sensor, 1 solve) | 69 (176) Dg; 52 Dg/fb | -0.38 / 0.66 | 1.32 1.24 / 1.29 1.25 / 1.22 1.10 | 0.01 -0.11; +0.509 +0.492 +0.627 | 0.89 1.64 0.34 / 1.22 -0.58 0.36 | 0.079 / 6e-17 | <1e-14 |
| rayflux (1 solve, linear-S rays) | 139 (176) Dg; 70 Dg/fb | -8.7e4 / 1.1e5 | 0.89 1.01 / 0.89 0.99 / 1.75 1.61 | 2.01 2.00; +0.001 +0.000 +0.000 | 0.43 0.66 0.63 / 2.05 1.59 0.72 | 0.088 / 2e-14 | <1e-14 |

Sources: `pos.out`, `pos_new.out`, `pos_rf.out` (breakdown: `python3 ana_pos.py <file>`),
`thin_0..3.out`, `thin_<variant>_<chi>.out`, `thick.out`, `thick_x.out`, `quick_new.out`,
`quick_rf.out`, `torder.out`, `torder2.out`, `oddeven.out`, `cons.out`, `rwave.out`, `picard.out`,
`front_alpha.py` (prints its table, about 1 min). Odd-even: the cfl 0.5 decay (0.1 -> 0.079 in 20
steps) is the physical compact-diffusion rate of a 2dx mode (exp(-0.13) for chi dx 100); A2 doubles it.

Other measured lines:
- Thick limit, extra alpha designs (`thick_x.out`): A1^2/1st L1 orders -0.99 -1.15 (kappa_num/kappa - 1
  +0.009 +0.017 +0.035, grows on refinement); A1^2/minmod 2.12 1.57; A1^2/mc 2.01 2.04.
- Time order, thin wave N 128 chi dx(64) 0.1, c dt/dx 4..0.25 (`torder.out`): SDIRK2 2.01 2.01 2.00 for
  npic 2, 3, 5 and for 1st order; backward Euler 1.12 1.07 1.04. Thick Gaussian N 128, c dt/dx 384..24
  (`torder2.out`): SDIRK2 + MC 1.99 1.79 0.84 (the last differences, 5e-7, sit at the Picard (npic 3)
  floor), BE 0.95 0.97 0.99.
- Radiating wave (`rwave.out`, `quick_new.out`; D fixed 1/3, |lambda_num - lambda| at N 32/64/128):
  ka 1: A1/1st 3.5e-1 1.8e-1 8.9e-2 (order 1.0); A1/mc 2.3e-2 5.9e-3 1.5e-3 (1.99 1.99); A3/vl 2.4e-2
  6.0e-3 1.5e-3 (2.02 2.00); S1 5.3e-2 1.1e-2 2.9e-3 (2.25 1.94). ka 100 (lambda = -0.10538): A1/mc
  1.1e-4 1.2e-4 5.1e-5; A3/vl 6.1e-4 4.2e-5 2.2e-5; A1/1st 8.6e-2 7.2e-2 4.9e-2; S1 5.9e-2 4.5e-2 3.8e-2.
  With ray closure the ka 1 mode converges to a different lambda (-0.913 + 3.117i), the transport
  (not Eddington) dispersion; no analytic transport reference was built. rayflux's Dfix lines are
  invalid (Dfix skips the rays, rayflux then read stale ray fluxes); its ray line at ka 100: 1.1e-3
  4.0e-4 1.9e-4.
- Conservation: sum E, sum F (vacuum) and sum (E + rho_cv T) drift <= 1.2e-14 relative for every
  variant run (all fluxes are face fluxes; stage 2 uses f(Y1) = (Y1 - U^n)/(gamma dt)).

### Realizability: what fails and why

No alpha candidate is realizable over the sweep; the violations come from four separate mechanisms,
identified by switching options one at a time:

1. Closure inconsistent with the state (dominant; pulses, counter pulses, slab reflection). The steady
   formal solution with an isotropic source cannot describe radiation that is streaming in the moment
   state but not tied to a boundary: a square pulse with F = cE in chi dx 0.1 gets D = 0.28-0.29 from
   the rays (D < f^2 = 1), and in the slab the rays put the reflected diffuse light into the whole vacuum
   region instantly (D drops 0.979 -> 0.69 while the moments still hold the pure beam, f^2 = 0.979).
   With D < f^2 the characteristic variable E - F/(c sqrt D) is negative, so the moment PDE itself
   produces E < 0 / abs(F) > cE; even alpha = 1 fails (pulse chi dx 0.1 cfl 0.5: -0.19 / 0.31). The floor
   D >= f^2 (`/Dg`, from the state's own flux, not from tau) removes these at c dt/dx 0.5 (alpha = 1:
   0 violations in pulse and slab), but it is not a rays-only closure.
2. SDIRK2 is not positivity preserving for c dt/dx > ~1: the stage-2 right-hand side
   U^n + (1 - gamma) dt f(Y1) = -1.414 U^n + 2.414 Y1 has a negative weight (no 2nd-order implicit RK
   is unconditionally positive). The BE fallback (`/fb`) cuts the worst violation from -0.31/0.69 to
   -0.036/0.146.
3. Closure Picard does not converge at c dt/dx ~ 4 (`picard.out` (a)): in the slab, even 1st order, the
   per-stage change stalls at ~1e-3 for 60 passes in 25 of 30 stages (2 of 16 at cfl 100; all converge
   at cfl 0.5 and 1000 within 4-18 passes). The thin smooth wave converges in 3-12 passes at every cfl
   for every limiter. The stall is the lagged D (rays + floor) switching, not the limiter (1st order
   stalls too).
4. Limiter (reconstruction) itself: slab at c dt/dx 0.5 with npic 1 / 3 / 10 / converged
   (`picard.out` (b), chi dx 10): 1st order -1e-97 (realizable) for all npic; minmod +1.9e-2 / 4.4e-2 /
   4.5e-2 / 4.5e-2; MC +8.1e-2 / 1.07e-1 / 1.07e-1 / 1.07e-1; van Leer +5.8e-2 / 8.2e-2 / 8.2e-2 /
   8.2e-2. Converged = npic 3 = npic 10: the slab |F| > cE at cfl 0.5 is NOT from lagged or unconverged
   limiting; the converged 2nd-order (unlimited-in-time, linear-in-the-iterate) reconstruction in the
   dissipation term is not realizability preserving at the slab front. At cfl 1000 all limiters are
   realizable to <= 1e-6 (chi dx >= 10).
5. Beams: A1 alone (no s_ray) gives abs(F) - cE up to 4.7e-3 on beams into absorbing media; A3 (alpha
   >= s_ray) removes every beam violation (A3/1st/Dg: beam 0, slab 5, all slab-scattering at cfl >= 4).
   RC weight 1 - theta changes nothing (86 vs 87); RC weight 1 (unweighted) is unstable in thin media
   (113 violations, worst -3e4).

Alpha at the steepest face (`front_alpha.py`, step 1 / step 4): A3 0.95-1.00 at every front and cfl
(slab at cfl 100 step 1: 0.10). rayflux 1.00 at the vacuum beam front but 0.00 at the pulse front in
chi dx 0.1 at cfl 0.5 and 4 (steady rays see no streaming there), and 0.00 / 0.02 at the slab front at
cfl 100 / 1000: rayflux drops the dissipation exactly where the moment state is ahead of or behind the
steady rays, and its sweep blows up (beam into chi dx 0.01 scattering at cfl 1000: -8.7e4 / 1.1e5).
S1: 1.00 at fronts at step 1 but 0.11-0.37 after 4 steps at cfl 100-1000 in the pulse (fronts smeared
by large steps look smooth to the sensor).

## Recommendation

alpha_f = max(A1(tau_f), s_ray_f) (A3) with RC weight 1 - alpha_f, SDIRK2, van Leer or MC
reconstruction of the dissipation-term states with 3 Picard passes per stage, the realizability floor
D >= f^2 and the BE step fallback.

- Order: it is the only family that is 2nd order in all three regimes measured here: thin wave 1.95-2.02
  (vl and MC), thick diffusion 1.95-2.21 with kappa_num/kappa -> 0 (vl +0.029 +0.007 +0.002; MC +0.003
  +0.000 +0.000), time 2.00 (thin), radiating wave 2.00. minmod is NOT enough in the thick regime
  (1.09, kappa_num/kappa halves per refinement) because alpha ~ 1/tau turns the O(dx^2) one-sided slope
  error into an O(dx) diffusion error; MC (or vl) is needed where alpha ~ 1/tau. MC is the most accurate
  thick choice; vl is AthenaK's PLM and nearly as good.
- Why not the linear alternatives: S1 is one solve per stage but measured 1.1-1.3 in the thin wave,
  0 in the thick limit (kappa_num/kappa +0.5, non-convergent: psi ~ 1 at the Gaussian peak gives
  alpha = A1 there) and 0.3-1.6 in time (closure and psi lagged to the stage start). rayflux gets the
  thick limit right (2.00, kappa_num -> 0, but only with linear-source rays), is 0.9-1.0 in the thin
  wave and 0.4-0.7 in time, and is unstable in transients (alpha -> 0 where the steady rays' flux
  disagrees with the moment state). Between vl + Picard and rayflux: vl + Picard.
- Cost: 3 Picard passes = 3 linear solves + 3 ray sweeps per stage (S1/rayflux: 1). In the thin wave
  the passes converge (change < 1e-12 in 3-12 passes at every cfl); 2 passes already give time order
  2.01 there.

Failures not fixed:
- Strict realizability (test 1) is not achieved by any variant: 42-60 violations remain with floor +
  fallback (worst -0.036 / 0.146), all in pulse / counter-pulse states (F = cE not tied to a boundary,
  where steady rays give the wrong D), the scattering slab at c dt/dx 4 (closure Picard stall), and
  thick absorbing pulses at the 1e-11 level. A realizable 2nd-order scheme would need either a closure
  that knows the moment state's streaming (time-dependent rays or the D >= f^2 floor, already used), a
  converged nonlinear closure solve at c dt/dx ~ 1-10, and a positivity limiter on the high-order part
  (the converged limited scheme violates at cfl 0.5 in the slab).
- The D >= f^2 floor uses the state's flux, so the closure is not rays-only; without it the pulse and
  slab tests cannot be realizable for any alpha (mechanism 1).
- Time order on the thick Gaussian flattens to 0.84 at the smallest steps (Picard floor at npic 3).
- No analytic transport (ray-closure) dispersion relation was built; test 5 validates the scheme with
  D = 1/3 only.

## Per-face positivity limiter (a posteriori, MOOD-like), closure rays-only (no D floor)

Implemented in `scheme.py` (`mood_stage`, `step` with `opt['mood']`), driver `mood.py`
(`pos [--all]`, `thin <variant> <chi>`, `quick`). Finished 16:55 PST 11-14; sweeps ran on vipa1235
(32 CPUs, allocation 12796548), the order tests on the login node.

Algorithm per SDIRK2 stage:
1. Solve the stage as before (A3, van Leer or MC, 3 Picard passes, rays-only D).
2. Mark cells with E < -1e-12 E_max or |F| - cE > 1e-12 c E_max; both faces of each marked cell
   (`mood_nb` = 1 adds one more face each side) become low-order faces: first-order states,
   alpha = 1, speed c, no Rhie-Chow. D and C1 of the cells are not changed.
3. Re-solve the stage with the marks; up to 3 rounds; marks accumulate over the step.
4. If a stage is still not realizable: redo the step with backward Euler and the same marking loop.
5. (`--all`, `mood_all`) If BE is still not realizable: the step in BE with every face low order.

Low-order flux on a marked face, two versions:
- `lf`: Lax-Friedrichs with a = c: F_E = (F_L+F_R)/2 - (c/2)(E_R-E_L),
  P = c^2 (D_L E_L + D_R E_R)/2 - (c/2)(F_R-F_L).
- `char`: upwinding of the characteristic variables E +- F/c with speeds +-c:
  F_E as above, P = c^2 (E_L + E_R)/2 - (c/2)(F_R-F_L), i.e. the face pressure average takes the
  free-streaming factor 1 on that face. Cell D, C1 and every unmarked face keep the rays' D.

Is the low-order flux with a = c realizable when the rays give D < f^2 (cause 1)? No, for any speed:
the LF convexity argument needs U +- Phi/a in the realizable set; for f = -1 this requires D >= 1
whatever a is (|-1 + (c/a) D| <= 1 - c/a only for D >= 1). Measured with every face low order,
backward Euler: pulse in chi dx 0.1 at c dt/dx 0.5 gives minE -0.121, max(|F|-cE) 0.233 (a = c) and
-0.063 / 0.155 (a = 2c); explicit LF at c dt/dx 0.2-0.4: -0.135 / 0.256 (a = c), -0.071 / 0.169
(a = 2c). The smallest measure that is realizable is the characteristic upwinding (`char`): every face
`char` in BE is realizable in all of pulse chi dx 0.1, counter chi dx 0.01, slab chi dx 10 and 1 at
c dt/dx 0.5, 4 and 100 (worst minE 0, max(|F|-cE) <= 0). It does change the pressure flux on the
marked faces (face factor 1 instead of the cell D), so it is not "D untouched" on those faces; it is
the minimal change found. With SDIRK2 instead of BE it fails at c dt/dx 4 (-0.047 / 0.286): the
stage-2 negative weight again, which the BE step handles.

### Realizability sweep with the limiter (176 cases incl. c dt/dx 1000; `mood_pos.out`, `mood_pos_all.out`)

| variant | violations | worst minE / max(abs F - cE) | marked faces per step | re-solves per stage | BE steps | all-faces BE steps | BE still bad |
|---|---|---|---|---|---|---|---|
| vl, lf, steps 1-4 | 49 | -0.39 / 0.55 | 70.3 | 1.32 | 37.7 % | - | 1904 |
| vl, char, steps 1-4 | 44 | -0.14 / 0.69 | 40.7 | 0.94 | 28.7 % | - | 817 |
| mc, char, steps 1-4 | 41 | -0.14 / 0.69 | 42.3 | 0.94 | 29.1 % | - | 772 |
| vl, char, nb 1, steps 1-4 | 42 | -0.057 / 0.71 | 61.5 | 0.73 | 24.9 % | - | 534 |
| vl, char, steps 1-5 | 3 (round-off) | -6.8e-13 / 2.2e-12 | 24.5 | 0.73 | 24.4 % | 369 (5.6 %) | 0 |
| mc, char, steps 1-5 | 3 (round-off) | -5.4e-13 / 2.2e-12 | 26.4 | 0.73 | 24.9 % | 375 (5.6 %) | 0 |
| vl, char, nb 1, steps 1-5 | 2 (round-off) | -6.8e-13 / 3.0e-12 | 48.2 | 0.62 | 23.0 % | 300 (4.5 %) | 0 |
| 1st, char, steps 1-5 | 2 (round-off) | -2.7e-13 / 3.1e-12 | 21.1 | 0.38 | 17.1 % | 352 (5.3 %) | 0 |

(6644 steps in total per variant; the remaining "violations" are 1-3e-12 of the reference energy,
the linear-solver round-off just above the 1e-12 test tolerance, in vacuum / near-vacuum pulses.)
More marking rounds without step 5 do not help: 10 rounds, slab chi dx 10 c dt/dx 0.5 still
max(|F|-cE) 0.117 with 2362 marks in 120 steps; the rays' D < f^2 fills the whole vacuum region, so
local marking chases the violation across it. The BE and all-faces rates above are for this
deliberately hostile sweep (square pulses with F = cE inside media, beams hitting slabs, c dt/dx up
to 1000).

### Orders and conservation with the limiter on (`mood_thin_*.out`, `mood_quick.out`)

The limiter never fires in the smooth tests (marked faces 0, re-solves 0, BE 0 in every run), and
the results are identical to the unlimited scheme:
- thin wave, c dt/dx 0.5, chi dx(64) 0.01 / 0.1 / 1: vl 2.01 2.01 / 1.95 2.00 / 2.01 2.02; MC 1.98
  2.01 / 1.97 2.01 / 2.01 2.02;
- thick Gaussian: vl 1.95 2.04, kappa_num/kappa - 1 +0.029 +0.007 +0.002; MC 2.21 2.20, +0.003 +0.000
  +0.000 (both dt ~ dx and dt ~ dx^2);
- time order (thin, c dt/dx 4..0.25): 2.01 2.01 2.00 (vl and MC); thick Gaussian 1.99 1.79 0.85;
- radiating wave (D = 1/3), ka 1: vl 2.02 2.00, MC 1.99 1.99;
- odd-even: 0.1 -> 0.079 (cfl 0.5), 3e-17 (cfl 100) in 20 steps;
- conservation with the limiter active (pulse tests: 11352 marked faces, 156 BE steps): sum E,
  sum F, sum (E + rho_cv T) drift <= 4.3e-15.

### Recommendation for the AthenaK port

A3 (alpha = max(A(tau), s_ray), RC weight 1 - alpha) + van Leer (AthenaK PLM; MC if the thick-limit
accuracy matters more) + 3 Picard passes + SDIRK2, rays-only closure, plus the a posteriori limiter
with the characteristic low-order face flux (`char`):
- per stage: solve, mark the faces of non-realizable cells, re-solve with those faces low order,
  up to 3 rounds;
- if still not realizable: the step in backward Euler with the same marking;
- if still not realizable: the step in backward Euler with every face low order (always realizable
  in the sweep).
The low-order face flux needs only a per-face flag in the flux kernel (alpha = 1, speed c, no RC,
first-order states, pressure average with factor 1); marking is a cell check plus a face mask, and
the re-solve reuses the same matrix assembly. The second-order results are untouched wherever nothing
is marked. `lf` (a = c with the cell D) is not sufficient: it cannot be realizable while the rays give
D < f^2.

## Staggered variant: E at cells, F at faces, F equation at every face (`stag.py`, `stag_tests.py`)

Finished 19:20 PST 11-14. The realizability sweep ran on vipa1235 (32 CPUs); the other tests ran on
the login node.

Scheme: `stag.py` docstring. E flux at faces Phi_f = F_f - alpha_f (a_f/2)[E]_f; F "flux" at cells
Psi_i = c^2 D_i E_i - alpha_i (a_i/2)[F]_i; van Leer jumps (lagged weights), alpha = max(A(tau),
s_ray) at faces and at cells, a_f = c sqrt(max D), a_i = c sqrt(D_i), no Rhie-Chow. Same rays (on
E and F interpolated to cells), same SDIRK2 + 3 Picard passes. The matrix is recovered from colored
probes of the affine operator (checked equal to the dense Jacobian to 0.0 for periodic and inflow).
Inflow: n+1 face unknowns; the left ghost cell is the characteristic split.
Positivity limiter: realizability is tested co-located, at cells (E_i, F interpolated to the cell)
and at faces (F_f vs E at the face, mean of the cells or the upwind cell). A bad cell marks itself and
its two faces, a bad face marks itself and its two cells. Marked locations get alpha = 1, first-order
jumps and speed c, with factor 1 in the pressure of a marked cell ('char'). Then the same chain as
the cell scheme: 3 rounds, BE with marking, BE with everything low order.

| test | staggered (vl) | cell (A3 + vl, same limiter chain) |
|---|---|---|
| thin wave chi dx(64) 0.01: L1(E)+L1(F), orders | 4.31e-4 9.82e-5 2.25e-5; 2.13 2.13 | 8.21e-4 2.04e-4 5.07e-5; 2.01 2.01 |
| thin wave chi dx(64) 0.1 | 3.84e-5 9.08e-6 2.21e-6; 2.08 2.04 | 9.69e-5 2.51e-5 6.27e-6; 1.95 2.00 |
| thin wave chi dx(64) 1 | 8.48e-10 2.11e-10 5.19e-11; 2.01 2.02 | 8.42e-10 2.10e-10 5.15e-11; 2.01 2.02 |
| thick, L1 orders; kappa_num/kappa - 1 | 1.94 2.02; +0.029 +0.007 +0.002 | 1.95 2.04; +0.029 +0.007 +0.002 |
| thick, first order (alpha = A(tau), no recon) | +0.866 (same fixed fraction) | +0.866 |
| time order thin (c dt/dx 4..0.25) | 2.02 2.01 2.00 | 2.01 2.01 2.00 |
| time order thick Gaussian (c dt/dx 384..24) | 2.02 2.01 2.00 | 1.99 1.79 0.85 |
| radiating wave D = 1/3, ka 1: abs(err) N 32/64/128 | 9.6e-3 1.8e-3 4.1e-4 (2.40 2.16) | 2.4e-2 6.0e-3 1.5e-3 (2.02 2.00) |
| radiating wave ka 100 | 1.4e-3 2.7e-4 3.9e-5 (2.36 2.79) | 6.1e-4 4.2e-5 2.2e-5 |
| odd-even (cfl 0.5 / 100, 20 steps) | 0.078 / 3e-17 | 0.079 / 3e-17 |
| conservation (sum E; sum F in vacuum) | 1.4e-15; 1.0e-14 | 4.3e-15; 2.3e-15 |
| realizability sweep, 176 cases | 88 violations (cell/face-mean tests); 62 on minE + upwind-face test | 3, all round-off (1-3e-12) |
| worst minE / cell / face | -0.021 / 0.235 / 0.380 | -6.8e-13 / 2.2e-12 |
| limiter work | 214 marked/step, 2.16 re-solves/stage, BE 64 % of steps, all-low 4190, still bad 3883 | 24.5 marked/step, 0.73 re-solves/stage, BE 24 %, all-low 369, still bad 0 |
| beam into vacuum, front 0.9 -> 0.1 width at t 0.4 (cfl 0.5 / 4) | 29 / 14 cells | 8 / 29 cells |
| beam onto absorber chi dx 10: E in the first slab cell (cfl 0.5 / 4) | 2.6e-3 / 1.0e-2 | 5.9e-7 / 9.2e-5 |

(The limiter never fired in the staggered smooth tests: 0 marked faces, 0 re-solves.)

Why staggered realizability fails:
- At a sharp front the staggered state is realizable as transport but not under any co-located test:
  F on the face behind the front cell is the beam value, E in the front cell is still filling, so
  (F_L + F_R)/2 > cE_i there. In the beam into vacuum this flags the front cell at every step
  (0.0165 at step 1, 3e-4 after 120 steps), even with every face and cell low order. The upwind-face
  test passes for that beam, but the mean-face and cell tests do not.
- The low-order staggered flux is not a characteristic upwind scheme: E and F live at different
  points, so E +- F/c are not upwinded. Even with everything low order the step stays non-realizable
  in 3883 of 6644 steps. In the slab the limiter makes the front cells low order at every step,
  which is why the cfl 0.5 front is 29 cells wide (cell scheme: 8) and the absorber edge leaks
  2.6e-3 (cell scheme: 5.9e-7).
- Thick absorbing slabs fail the cell test (+7.7e-3) because interpolating F across a cell inside an
  exponential decay of e^-10 per cell overshoots cE_i.

### Recommendation (staggered vs cell)

Keep the cell scheme (A3 + van Leer + Picard + SDIRK2 + rays + per-face limiter with 'char').
- Accuracy per cell: in smooth fields the staggered scheme is better. Thin-wave error is 2-2.5x
  lower, radiating-wave error 2.5-4x lower at ka 1, and its thick time order is clean (2.00 vs the
  cell scheme's 0.85 floor at small steps). The thick limit is the same: both give second order and
  kappa_num -> 0 with van Leer. Staggering needs no Rhie-Chow.
- Positivity: the staggered scheme cannot be made realizable with the tools tested here. There is no
  unambiguous co-located test, and its low-order fallback is not realizable, so 88 of 176 cases
  violate with the limiter, against round-off only for the cell scheme. Its limiter work is about
  9x that of the cell scheme (BE in 64 % of steps), and sharp fronts and absorber edges come out
  worse in accuracy, not better.
- Simplicity: the staggered operator is simpler (no RC, compact diffusion for free). The positivity
  machinery is harder: two marking sets, interpolation choices, and a different boundary F unknown.
The staggered layout would only pay off if a realizable low-order staggered flux and a realizability
criterion that accepts transported fronts were found. Neither exists in this prototype.
