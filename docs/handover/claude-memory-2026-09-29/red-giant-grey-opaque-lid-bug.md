---
name: red-giant-grey-opaque-lid-bug
description: "LATENT BUG found 2026-09-09 08:20: on the grey path hot_jupiter_param.grav (and ap, Rgas...) were set only under rt_ck, so grav = 0 and the top ghost's dtau = kappa p/g_eff = +inf -- an OPAQUE LID above the domain in EVERY grey run to date, prod11 included. Fixing it doubles L_rad,out at t=0 (0.5 -> 0.9-1.0 L) but the thin atmosphere then cools 3364 -> 2200 K and contracts (the old top-cooling runaway?). Being isolated on the plain 1.1 R star (R8_lid_open / R8_lid_wall)."
metadata:
  type: project
---

Found by the corona agent when kappa -> 0 in the corona turned the lid's inf into 0/0 = NaN.
`src/pgen/red_giant.cpp` ~883-901 now sets hot_jupiter_param.{grav, ap, Rgas, ...} on the
grey path too.  Consequences measured in R7_hot_open (8x8, 1.5 R, corona at 6e5 K):
- L_rad,out/L at t = 0: 0.87 -> 0.97 (R5 with the lid: ~0.5).  So a large part of the
  "0.31-0.6 L resolution deficit" ([[red-giant-flux-deficit-is-spinup]], prod11's 0.31)
  may be the LID back-warming the atmosphere and suppressing the emergent flux.
- The thin atmosphere (tau < 0.1) then cools and contracts: T(i=300) 3367 -> 2256 K by
  2e5 s, rho there down 33x (hydrostatic shrinkage of a cooler skin); the hot corona,
  itself inert and hydrostatic at t=0 (H = 1.56 r, join pressure-matched to 1e-5, wall
  clean, deep quiet, Mdot out 4.5e-21 Msun/yr), drained into the vacated space at its own
  sound speed and the run died at 3.09e5 in the photosphere region (i 239-311).
- Physics or runaway?  A cooler skin is expected if < L emerges (T_skin ~ L^(1/4)), but
  L_out was RISING toward L while the skin fell to 2200 K -- not RE.  Matches the old
  [[red-giant-top-cooling-runaway]] that the accidental lid has been masking since.
- The corona needed dfloor 1e-26 and eos_logd_min -27 (rho 6e-25 pressure-matched), and
  the conduction/RT cutoff <hydro>/rad_kappa_rmax + rad_kappa_above (kappa above it; 0)
  -- both implemented (conduction.hpp/.cpp, two_stream_rt.hpp) -- with the cutoff at
  3.885e12 (the join face 3.8871e12; 3.90e12 would leave two corona cells radiative and
  let their K ~ 1/kappa set the dt).

**prod11 runs on the pre-fix binary: it HAS the lid.**  Stable, but its 0.31 L and its
photospheric T are partly artefacts of the lid.  Do not restart production until R8
settles whether the lid-free atmosphere equilibrates (then production should be restarted
fresh with the fix -- radiating ~L) or runs away (then the thin-cell RT needs work first).

## R8 RESULT (09:10): lid quantified, NO runaway, but lid-free is LESS stable at 8x8
- Top ghost old vs new: dtau_top +inf -> 1.8e-5, I_down/B 0.846 -> 2e-5; rt_top_re's source
  is I_up/2 so the lid halved the net top flux EXACTLY: L_out/L at t=0 0.4994 -> 0.9989.
- Thin cells track T_skin = 4000 (L_out/L)^(1/4) / 2^(1/4) up and down; d ln T/dt at i=300
  alternates and decelerates; the scheme's thin-limit equilibrium is T_ph/2^(1/4) exactly
  (emission/absorption 0.96-1.18 numerically).  The old "top-cooling runaway" is NOT here.
- Evolved L_out/L is 0.35-0.57 WITH OR WITHOUT the lid -> the deficit is sub-photospheric
  transport (resolution), not the top boundary.  Lid-free t=0 emission 1.0 L decays to
  0.37 by 1e5 s.
- Stability: R8_lid_wall (prod11's config + fix) died at 2.5e5, domain-wide NaN (the wall
  ghost holds the INITIAL 3364 K column while the lid-free top cools to ~2600 K -> the
  pinned ghost fights the interior; the pgen comments warn of this).  R8_lid_open died at
  1.307e6 on a DEEP conduction collapse at i~137 (r/R 0.6) after dt degraded 59 -> 3 s
  (being localized: real mid-envelope event or a table/floor artefact?).  R6_open (lid)
  completed 1.5e6.  => prod11 is stable BECAUSE the lid holds its top at the IC T.
- Next: R9_lid_open32 (6 nodes, 32x32, prod11 config + grav fix + outer_bc open, tlim 2e6):
  if stable and sane it replaces prod11 as production (fresh start; outer_bc cannot be
  switched on a restart).

## Localized (09:30): R8_lid_open died of MY opacity floor, not physics; R8_lid_wall of the wall
- R8_lid_open: cell (0,9,5,279), r/R 1.063 (the PHOTOSPHERE, one block corner), rho 1.75e-8,
  T 5.92e6 (7136 K one dump earlier -> an 800x jump between dumps), tau 65-281, w 0.77-1.0,
  logR -10 -> below the window -> got the thin-gas floor kappa 1e-5 -> kappa_rad 3.6e29 cgs
  -> conduction dt 1.4e-7.  Shell i=120-160 inert and identical in R8/R6.  FIX applied:
  the logR floor now applies ONLY to the prepended rows (rho < 1e-14); the file's own
  below-window nodes keep their edge-filled values (physically: e-scattering/free-free
  for a hot shocked cell).  Rebuilt; R9 relaunched on the fixed binary.
- R8_lid_wall: its reflecting ghost held the IC (3364 K, 1.5e-11) against a top that had
  cooled to 2100 K and thinned 6x, all falling at -4e5 coherently -> domain-wide NaN in one
  cycle (dt = NaN from a bad cell).  The WALL is incompatible with a lid-free top.
- R9_lid_open32 (job resubmitted ~09:35 after a launch race on the binary path): 32x32,
  prod11 config + grav fix + outer_bc open + floor fix, tlim 2e6.  If stable and sane it
  replaces prod11 as production (fresh start).

## 10:30 R9 post-mortem: the open ghost is CLEAN; lid removal is the variable
Top three rows at 1e5: hydrostatic continuation uniform to 0.1 % over all 6144 columns,
no failed inversion / sentinel / out-of-table state, zero outflow; extremes in plain
interior columns; prod11's lidded top is MORE extreme at 2e5 (v_r -1.4..-1.9e5, T 2367,
rho 8e-12) and survives.  R9 died like the 8x8 lid-free WALLED star (domain-wide NaN,
1.99e5 vs 2.5e5): open bought nothing.  dt flat to the last print, then one step.
Mechanism still open: a one-step all-columns NaN in the RT domain (i >= icut ~ 200) of
every column -> the RT step hitting a shared top state.  Prime suspect: MY equilibrium
closed form's `absn <= 0 -> deq = -e` branch (relax to zero energy), reachable only with
no downward intensity.  4-way bisection at 32x32 (tlim 5e5, ~10 min each, jobs 194562-5):
T1_wall, T2_semilin (old linearized step), T3_nonewton, T4_topre_off.  Monitor reports in
one table.  Ghost guard (corona agent) is hardening only, not the cause.
10:45 user: "maybe the semi-implicit approach is not correct after all -- try a full explicit
test."  New problem/rt_explicit (de = src*dt, no relaxation/Newton; stability limit = local
radiative time ~1e3 s at tau 1-10 vs dt 30 s) being added by the corona agent together with
the open-ghost guard, the RT top ghost from the ACTIVE cell, and a post-Newton positivity
guard (the Newton loop could exit with ei + de <= 0 -- a real one-step NaN path).
T5_explicit (32x32, R9 config + rt_explicit, tlim 5e5) is staged in red_giant/T5_explicit,
to be launched on that binary alongside T1-T4.
