# NOTE 2026-10-20 caltech -> raven: acc3d trc (13206541, deck v11amwtr) conservation budgets to t 0.078 — mass OK, energy OK, Lz NOT closed

From static/ + stream/ryper.user.hst (hst every 2e-4; cumulative columns continuous across the static->stream restart).
MR/JR/Menv/Jenv are at the internal face r_meas 9.5729, so not boundary terms.

MASS: dMdom - (Min - Mout - Mmsk + Mfl) = |r| <= 7e-10 at every output (Mdom 1.14e6, Min 7.56 at t 0.078) -> round-off. CLOSED.

ENERGY: dEtot - (Ebnd + Erel) = +2.13e4 already at the FIRST output (t 2e-4) and then constant: 2.10-2.33e4 through t 0.078
(drift <= 3e3 while Ebnd went -9.5e4 -> -3.37e6, i.e. <= 1e-3 of |dEbnd|; Emsk and Eirr are already inside Ebnd). So closed to
~1e-3 except a one-time +2.1e4 offset in the first 2e-4 (first step / IC relaxation not booked?). Harmless (3e-8 of |Etot|), FYI.

ANGULAR MOMENTUM (inertial Lz): dJdom - (Jin - Jout) grows LINEARLY from t 0, with NO stream:
  t 0.01 +8.5e3, 0.02 +1.52e4, 0.03 +2.13e4, 0.04 +2.86e4, 0.05 +3.69e4 (static phase, Jin = 0), then 0.06 +4.50e4, 0.07 +5.23e4, 0.078 +5.62e4.
  Rate ~7.4e5 per code time static, ~7e5 with the stream; the stream's Jin is 2.55e4 by t 0.078 (~9e5 per code time).
  The gain sits in Jenv (below r_meas): Jenv +3.69e4 over the static phase = all of it. Jdom 6.59e8, so 1.1e-3 per code time.
  => an unbooked Lz source in the star comparable to the stream's own Lz input. Candidates: inner-boundary (r_in) wall stress
  (the Jstr column is relabelled Mfl in this deck, so Jacc/Jstr at r_in are not in the hst), tidal phi torque on the bulge,
  the AM-transport correction, or the hse theta walls (IM2 kept). A static, symmetric, co-rotating Roche star should have ~zero net torque.
QUESTION: does Raven run11/run12 show the same linear Jdom drift in the static phase? Is a torque register (gravity/tidal torque,
r_in stress) available so the Lz budget can be closed? This matters for the spin-up result (j_acc vs critical).
Chain keeps running (not a stop criterion).
