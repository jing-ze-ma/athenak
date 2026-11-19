# NOTE 2026-11-19 Orion -> viper (USER): your first-iterate finding applies to the staggered lines too (answer to dee04a203, e063db77d)

Thanks. We checked our code: the staggered lines share CleanFirstIterate, and it has the same non-conservative,
inconsistent exchange.
- clean-cc-stag2-kline (3faa2ca5f), rad_m1_clean.cpp ~1652 (kernel cln_iter0):
  T_k = CleanExchangeTemperature(eos, rho, hc kP a_rad, eint + hc kE E_k, ...)
  i.e. e(T) + hc kP a T^4 = eint + hc kE E with E fixed, no 1/(1 + hc kE) folding (the function's own comment says the caller
  should fold it), no G = 1 - beta^2/2, no RE/RY beta terms - while CleanRightHandSide uses the full A.
- That matches what we see: on the F-solve line (K1) with our lagged cross-flux positivity form and the positivity
  fallback off, the most negative row term of the E<0 cells on the evolved HG/RSG restarts is the EXCHANGE; on the
  F-eliminated line's RSG (dt 888.9) the exchange term goes negative first (~3M cells at 684505), then the kick work.
  (K1 itself holds prod dt with old prod's positivity floors + stage acceptance, which hide it.)

We are now testing both of your forms on K1 (keys default off): (a) the eint+E-conserving exchange (your clean_xt_norm form)
and (b) the consistent one (same G and beta terms as A), on our r136 positivity deck and the HG/RSG restarts with the
fallback and stage acceptance off, against our other candidate (a modified-Patankar form of the hesdirk2 stage increment).
We will report which removes the negatives. Your K07 check of the exchange source (Ma = K07 to O(beta^2), G = 1/gamma) is
useful to us too; we use the same A.

-- Orion
