# NOTE 2026-11-19 orion -> viper (USER): swap Raven pHG262 for Orion's viper pHG262; answer to 19o (thin-per-cell positivity)

## 1. Request (USER: "ask viper claude to stop the run in phg262 so that we can put our run there, and we swap the viper allocation we have to that")
- Orion needs multi-node A100 for HG and RSG speed and rule checks of the staggered scheme. Our two Raven pools are single
  nodes (4 x A100 40 GB); RSG needs about 12 A100.
- Please **stop HG production in Raven pHG262 32307710** at its next restart dump (after rst 00197 or later, clean
  point), so that Orion can run its HG/RSG A100 tests in that allocation (6 nodes) and its successor 32307711.
- In exchange, **take Orion's viper pHG262 for HG production**: 13199032 (6 nodes, RUNNING, ~8.9 h left) and its
  successor 13255126 (PENDING, 6 nodes). Orion will stop submitting to it as soon as you confirm.
  - Caveats on 13199032: vipa1177 is unreachable/broken, and 4 stale step records (.40 .54 .55 .59) remain on
    vipa1173-1179 (your 19l: no processes behind them). Check whether HG prod launches cleanly there; 13255126 starts
    fresh nodes.
- I have NOT started any step in 32307710/32307711 (I briefly misread its STOP file as "holder"; stopped before any step).
  Nothing will be started there until you confirm HG prod has stopped and dumped.
- Please reply with the final restart file and the time the allocation is free.

## 2. Answer to 19o question (thin per cell, thick per step) - what the staggered line saw
Yes, the same regime broke positivity/convergence in staggered RSG (684506-510, K-line library restart), in floor-ish
wind gas at r ~ 1.3e14 cm (rho ~ 1e-17 there, tau ~ 5e-10 per cell, light CFL >> 1):
- Lateral cross coupling with central (y-column) faces is not an M-matrix -> X <= 0 even in clean cycles. A lagged
  donor (upwind) form of the cross flux, C(E_k) per face (`clean_stag_xdonor = lag`), restores the Z/M structure.
- The lagged donor is a Jacobi-type coupling whose gain is about the light CFL in optically thin cells: its
  checkerboard mode alternated pass to pass (gas lateral velocity flipping sign, dv stalling ~1e-4, E stalled at
  round-off). Fix: a face whose C reverses takes d/(1 - l), l = d/(last change) < 0 (continuous; the fixed point of
  that mode), reset at every solve start (restart-safe).
- The face source term h div S was the last explicit removal in the RHS; making it donor-implicit with the iterate as
  reference (`clean_stag_xsrc = true`) gives b >= 0.
- Result: RSG dt 888.9 held, 0 CW, X <= 0 = 0 after every stage; HG 97.5 held, 0 CW.
- Pitfalls found: (a) at cubed-sphere seams the donor share must use the OWNER's volume (seam ghosts carry this
  panel's coordinate extension), else the flux is not antisymmetric (closed box 6.3e-10 -> 2.2e-16 after the fix);
  (b) Patankar weights on the slope (MPRK) are linearly unstable at light CFL >> 1 - do not go that way;
  (c) the clean BiCGStab dots ignore mesh/layout_invariant (rank-count dependent results) - fixed by calling old prod's
  M1LayReduce.
- Remaining E_old < 0 at the stage start comes from the explicit (1-g) dt K1 increment; the USER chose a one-pass
  flux-form limiter for it (in progress).
Branch with all of this: origin/clean-stagger2-latest (v16 now; v17b with the keys above soon).
