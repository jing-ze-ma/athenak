---
name: pcr-timing-result
description: PCR (parallel cyclic reduction) column solver timings 09-16, target 1.2x standard NOT met; best 1.73x (nseg32), 1.47x with maxit2/tol1e-6
metadata:
  type: project
---

Bench bench/bstar_fecz/m3pcr (branch rg-box-pcr, b56c30b1). Wide cost run 336x64x64, 200 cycles, GPU, cpu-s: w_std 26.92, w_thomas 67.44 (2.50x), pcr16 47.6 (1.77x), pcr32 46.5 (1.73x, best), pcr64 53.4 (1.99x), pcr64 + maxit2/tol1e-6 39.4 (1.46x), +cvfreeze2 no gain. Gate: thomas (T) bitwise = ref (R). PCR vs thomas differs at 1e-4 in user.hst (not bitwise, expected). PCR nseg=1 gate (P1) never finished. arm.sh arms T2/T3/P/P2/S lack the prod COMMON switches (uform, kiter, sts_perplane, tr_tau) so their timings are not prod_m3-comparable.

**Why:** the user's bar for using PCR (mode 3 <= 1.2x standard) is not reached; the cost is the Newton passes and horizontal operator, not the tridiagonal solve alone.
**How to apply (superseded 09-16 by the user): PCR (nseg 32) GOES INTO the prod_m3 chain at the re-gate after the bottom-flux fix; binary = wt_m3acc build_hip_pcrfix (fix + PCR), gate arms G3p/G0p in final_gate, staged input prod_m3/fecz_rt.athinput.pcrfix.** Earlier rule: the decision (accept 1.5x with loosened tolerance for He/rg/dhj, or drop) is the user's. See [[bstar-column-pulsation]], [[implicit-transverse-raddiff-plan]].
