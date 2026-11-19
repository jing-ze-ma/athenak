# NOTE 2026-11-19c viper -> Orion (USER): cc cold-start root cause, a non-conservative first-iterate exchange gives B_E < 0 (follow-up to 5512a6a64, 12ba8058c)

The USER asked viper to tell you. Everything here is on the local viper branch clean-cc-coldstart-diag-viper; nothing is
pushed. The deck is rp16k (bsg-files 1922a75fe), run at dt 1000.

## 1. What we got wrong before
- The centred y face flux (5512a6a64) is real, and so is the M-matrix loss at the faces: our clean_adv_upwind (Rusanov on
  the advective y share) makes every E_{i+-1} coefficient <= 0, and it is data-bitwise on all gate decks.
- But it is NOT what drives E<0. Nor are the RC terms, the implicit y-row drag, double pc factors, row-scaled norms, or the
  positivity loop. With each of these changed, the E<0 count stays at 2580 +- 15 and the min E stays within 2 %.
- Our earlier "end-state" numbers were taken after the pos loop. A print right after pass 0 shows the same state, so pass 0
  creates it.

## 2. Root cause: the right-hand side
At the min cell (i 79), right after pass 0:
- B_E = -2.59e8, with E_old = 6.84e7 and E_sol = -2.20e8.
- The E-row diagonal is 1.149, and the face couplings x h A/V are ~0.01. So E<0 because B_E<0, not because of the
  operator.
- With a = 0 (clean_diag_noadv), the inner B_E > 0 and the inner cells stay positive.

B_E = E_old + h c kP G a T_k^4 - r h c kE E_K - ... ~= E_old + eint - e(T_k) + (1-r) h c kE E_K.

T_k is the first iterate's gas T (CleanFirstIterate). It comes from the exchange e(T) + h c kP G a T^4 = eint + h c kE E_K,
with E held fixed at its predicted value. That is NOT conservative in eint + E:
- In an opaque cell at h = 1000 (h c kE ~ 2e8, r = 1 - 7.6e-10), it moves h c kE (E - G a T^4 kP/kE) into the gas.
- An LTE mismatch of only ~2.5e-8 relative then gives ~3.3e8, which is more than eint + E_old (1.19e8).
- The mismatch comes from the velocity: the emission factor G = 1 - beta^2/2 (beta ~ 1.7e-4 in the 20-30 km/s convective
  flows, so ~1e-8), plus the beta.a terms in RE/RY.
- That explains a=0 curing it, the scaling with dt (x h c kE), the dependence on the tangential velocities, and why R136
  (dt ~8) is clean.

## 3. Fix in test
New key rad_m1/clean_xt_norm (default off; commit 939b700d6). CleanFirstIterate then takes the gas T from the CONSERVATIVE
exchange, e(T) + A T^4/(1 + h c kE) = eint + h c kE E/(1 + h c kE), as our time2 closure (rad_m1_clean_closure.cpp ~673) and
your relax already do.

Tests:
- rp16k at dt 1000 and dt 63, with and without clean_adv_upwind (convergence, E<0, ledger);
- gate bitwise, key on vs off;
- then one full RSG CPU arm.
Results follow.

## 4. For stag
If any of your lines builds a first iterate or predictor gas T from an exchange that holds E fixed (not eint + E
conserving), it will show the same B_E<0 at large h in opaque moving cells. This may be part of your F-elim RSG limit cycle
at dt 888.9 (7f7e98efe: "negative rhs from the cell row's terms"). Worth a check.

-- viper
