# NOTE 2026-11-19d viper -> Orion (USER): the clean exchange source checked against Krumholz et al. 2007; first-iterate inconsistency (follow-up to dee04a203)

The USER asked viper to check our exchange equations against Krumholz, Klein, McKee & Bolstad 2007 (ApJ 667, 626;
astro-ph/0611003) and to tell you the result.
- We read the paper via ar5iv, Eqs. 6-52. The equation numbers were not checked against the PDF.
- We could not find Ma et al. 2025 online, so the check is on the form the code writes (rad_m1_clean_cell.cpp ~765,
  "Ma et al. 2025 eq. 56").

## 1. The code's form
h S_E = A - W:
- W = v.dm, the work of the radiation force; it goes to the gas KE.
- A goes to the gas eint:
  A/(hc) = (kP a T^4 - kE E) G + kE beta.(y + F'/c), with G = 1 - beta^2/2, y = F/c, F' = F - aE and aE = vE + v.P.
- Equivalently, A/(hc) = kP G aT^4 - kE E + RE E/(hc) + RY.y/(hc), with RE = hc kE (beta^2/2 - beta.a/c) and RY = 2 hc kE beta.

## 2. Comparison with Krumholz 2007
K07 Eqs. 16-17 are exact (Mihalas & Auer 2001); Eqs. 24-25 are the O(beta^2) expansion. The radiation obeys
dE/dt + div F = -c G^0, and the gas momentum gets +G.

The radiation source must equal -c G^0, so A = -c G^0 + v.G. Substituting Eqs. 24-25 with kE = kP (K07 Eq. 21):
A/c = kP(aT^4 - E)(1 - beta^2/2) + kP beta.(2F/c - aE/c)
This is the code's A, term for term. So Ma's form is K07 to O(beta^2), split into the heating (A) and the work (W).

Where G comes from: in the exact Eqs. 16-17 the emission enters G^0 as gamma kP aT^4 and G as -gamma kP aT^4 v/c. Its share
of A is therefore gamma (1 - beta^2) kP aT^4 = kP aT^4/gamma. So G = 1/gamma, exactly. K07's (1 + beta^2/2) on (aT^4 - E)
is the total gas-energy source; Ma's (1 - beta^2/2) is the heating alone.

Check: in a uniform moving equilibrium (E0 = aT^4, F0 = 0, E = E0(1 + 4 beta^2/3), F = 4/3 v E0), both forms give exactly
0.

Differences that are not errors:
- Ma keeps kE and kP separate, as in the exact Eq. 16.
- Ma keeps all the O(beta^2) terms, while K07's minimal set (Eqs. 40-42) drops some of them by regime.

## 3. What is inconsistent: our first iterate, not the equations
- CleanFirstIterate (rad_m1_clean.cpp ~1140) sets T_k from e(T) + hc kP aT^4 = eint + hc kE E, without G and without
  the RE/RY beta terms.
- CleanRightHandSide uses the full A.
- In opaque moving cells the two balances differ by O(beta^2, beta.a/c, beta.y/E), about 2-3e-8 relative.
- Since r = 1 - 1e-9, the T^4 terms of B_E do not cancel, and the mismatch times hc kE (about 2e8) gives the
  B_E = -2.6e8 of dee04a203.
- Two fixes are in test (local, not pushed):
  - (a) clean_xt_norm: the eint+E-conserving exchange, which divides the mismatch by 1 + hc kE;
  - (b) a "consistent" variant: the first-iterate exchange uses the same G and beta terms as A.

## 4. For stag
Please check that any exchange your lines solve for a predictor or first-iterate gas T uses the same G = 1/gamma and the
same beta terms as the stage's A. A bare kP aT^4 = kE E balance against a full-A RHS gives this B_E<0 at large h in
opaque moving gas.

-- viper
