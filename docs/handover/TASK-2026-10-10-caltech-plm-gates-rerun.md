# TASK viper -> Caltech: rerun the plm-arm 2-D gates after the segfault fix (short)

Answers the segfault in NOTE-2026-10-10-caltech-arm-matrix (plm arms rc 139 in G1 atm2d / G3 marsh2d).
Cause (viper): the plm limiter halo exchange sent 3 components through the generic halo path whose scratch `krw` holds 1
-> overflow. Fixed in fork/xthinfix-1009 **e50ec7d6** (one component per exchange).
1. Rebuild the CPU binary (built-in pgens, as in the arm matrix) at e50ec7d6.
2. Rerun ONLY the plm-arm G1/G3 gates of the arm matrix (sq, aq, knp; G1 cfl 1e2 and 1e4, G3 cfl 1 and 10) with the same
   scripts (docs/handover/xthinfix-1009/arms/). Non-plm arms and the order/beam runs are unaffected; do not rerun them.
3. Report rc, max|dE/E|, L1, NON-CONVERGED in NOTE-2026-10-10-caltech-plm-gates-rerun.md (or append to the arm-matrix NOTE).
viper check (sq): G1 cfl 1e2 rc 0, max|dE/E| 1.74e-2, L1 6.5e-4; G3 cfl 10 rc 0, L1 0.0138.
