# TASK 2026-10-16 caltech -> raven: accretor 3-D (+-20 deg) plan (FYI) + a short 3-D benchmark

**Test only, report back by NOTE** (docs/handover/NOTE-2026-10-1x-raven-accretor-3d.md on bsg-files-1009).
Branch `accretor-bigdom-1015` (github jing-ze-ma/athenak), head **08625ed7** (code 049c1f39). Design:
docs/dev/accretor_rhd_design.md **section 14** (14.1 says what is and is NOT 3-D yet: read it first).

## (a) FYI: the 3-D plan
The user wants a 3-D version of the bigdom accretor: theta range +-20 deg about the orbital plane instead of the
+-0.35 deg band. User's planning cost estimates: **~2,500 H200-h per orbit at full resolution +-20 deg on the old
domain, ~300 at half resolution on 4 GPUs, times ~1.3 for the big domain.** (Caltech's 200-cycle sanity run below
measured less, pre-impact: ~30 H200-h/orbit at half resolution on the big domain; your benchmark is the check.)

What 049c1f39 makes 3-D (switches on automatically when theta reaches > 2 deg from the plane; 2-D decks are
bitwise unchanged, gate 14.5):
- Stream nozzle Gaussian in z = r_f cos(theta) with sigma_z = **1.0054 Rsun** (problem/stream_sigma_z; the one-zone
  core at r 14.5 of stream_path_1015, same model as the 2-D deck); peak density normalised so the full 3-D stream
  carries **stream_mdot = 254.843 code = 1e-4 Msun/yr** (exact discrete sum over the 3-D window, |z| < 3 sigma_z).
  Its central slice is 18.01 in 2-D band units (the 2-D deck's 19.795 is 10 % denser than the one-zone core at
  14.5; alternatives in 14.2, a user decision).
- Donor mask = the exact 3-D Roche lobe (Phi(r, theta, phi) < Phi(L1), x > d_L1); donor_irr unchanged on it.
- Back sink and the Mnzf register restricted to the z window. Theta walls: **reflecting** (14.3).
- NOT yet 3-D: gravity, Coriolis, IC, WB and ambient are still the band's (theta-invariant; no vertical gravity).
  The 3-D deck is a COST AND PLUMBING benchmark, not physics.

## (b) Prerequisite
Your 2-D v3 test (TASK-2026-10-15-raven-accretor-bigdom.md, v3 1db9e3e7) passes first. If it has not, finish or
report that one before this.

## (c) Short 3-D benchmark (after (b))
Deck: **docs/handover/accretor-bigdom-1015/b3d_half.athinput** on the branch (same tables/column file as the 2-D
deck in that directory; `sed "s#@HDIR@#$PWD#" b3d_half.athinput > run.athinput`).
- Grid 420 x 114 x 1024 = **49,029,120 cells** (half resolution; theta 69.96..110.04 deg at dtheta = dphi =
  2 pi/1024), 4 MeshBlocks 420 x 114 x 256 = 4 GPUs (keep one block in x1). Stream on from t = 0, nlim 200;
  add `problem/dt_report=true time/ndiag=10`.
- Use Raven's own GPU/node rules (GPUs per node, partition/QoS, exclusions, GPU env exports, walltime, scratch).
  A few hundred cycles (nlim 200-400) is enough. Delete the rst files afterwards (7-11 GB each).
- Measure: s/cycle (steady part; the first ~120 cycles are a switch-on transient), dt and where it is set
  (dt_report), GPU memory, and **GPU-h per orbit** = (P_orb 0.473169 / dt) x s/cycle x nGPU / 3600 (say which GPU).
- Check: no FATAL / NON-CONVERGED / dt retries; Mnzf rate (hst col 25) = 254.843 within 0.5 % after the fill;
  nozzle z-Gaussian (hydro_w bin: first nozzle cell r ~14.51 at phi ~0.034 rad, rho vs z_f = 14.47214 cos(theta)
  against 0.2493 exp(-z_f^2/(2 x 1.0054^2)), floor beyond |z_f| 3.02); mask quiet (Mmsk col 9 small, smooth,
  lobe at floor at all theta); donor irradiation present at all theta (Eirr col 26 grows, m1 bin E = a T_d^4 in the
  lobe), Ltop sensible.

Caltech reference (049c1f39, 4 H200, 200 cycles, t 0.00282): no FATAL, 0 NON-CONVERGED, 0 retries, Picard mean 2.8;
dt 1.408e-5 flat, set by an envelope cell at rest (r 8.786); Mnzf 254.93 -> 254.89; Mmsk 7.4e-5; Ltop 1.44e5 ->
1.90e5; **0.80 s/cycle steady** (1.0-1.9 in the first 120 cycles); 24.6 GB per H200; ~30 H200-h/orbit pre-impact.

## (d) Scope
Test only: no long 3-D runs. Report back by NOTE (numbers above, plus anything that broke).
