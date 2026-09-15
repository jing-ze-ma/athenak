---
name: warm-started-newton-plan
description: NEXT after the chains are healthy (user 09-15): warm-started Newton for the mode-3 column (extrapolate the iterate from the previous stage/cycle; frozen Jacobian on RK stage 2) to cut the 2-3 passes toward 1 - the column solve is >50 % of the cycle after PCR + halo_every
metadata:
  type: project
---

Target: mode-3 column solve (RTCol3::Solve / RTCol3TeamSolve in wt_m3acc src/utils/two_stream_column_*.hpp), currently 2-3 Newton passes per stage at tol 1e-8; the passes are the cost (cvfreeze arm: assembly 1.5 %). Plan: (1) warm start = initial iterate extrapolated from the previous stage's converged column (store the converged b/T per column; at the next call start from it or a linear extrapolation in time) instead of the entry state; (2) frozen Jacobian on RK stage 2 (reuse the stage-1 factorization, one corrector pass); (3) count passes per call (it_mean is printed in the rt_col3 stats) - success = it_mean -> ~1-1.5 with the SAME converged answer (tol unchanged) -> results equal to roundoff, not bitwise; gate with analyze.py G3 numbers (Ftop 1.00000, floor v_rms, Mmax, P_nyq) + He Q1p2 fluxes. Expected 30-40 % off the column = ~20 % of the cycle. Radiation-to-hydro parity is NOT reachable (would need 10x); the only large-factor route is horizontal coarsening of the columns in the diffusive region (physics approximation, separate validation). Do NOT start before the B-star chain (11710002) passes its first turnovers and the He chain is swapped to 76b42ca8. See [[bstar-prod-cost-profile]], [[rkl1-halo-every]].
