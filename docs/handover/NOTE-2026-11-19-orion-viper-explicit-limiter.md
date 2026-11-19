# NOTE 2026-11-19 orion -> viper (USER "tell viper"): swap done on our side; one-pass limiter for the explicit increment

## 1. Swap (your 19p)
- Thanks. On Orion's side: gpuq daemon and the viper scheduler are stopped, the viper queue is disabled
  (vgpu/queue.txt.HANDED_TO_VIPER_1119). Nothing of ours will start in 13199032/13255126. Our last test step there had
  already ended; the 4 stale records (.40/.54/.55/.59) remain (no processes behind them).
- We now use Raven pHG262 32307710/11 only via --overlap steps.

## 2. Positivity of the explicit increment (USER decision: conservative, counted limiter for the explicit part;
## implicit stage solve stays positive by construction)
Staggered line, hesdirk2, explicit increment (1-g) dt K1. Findings (S3P8, branch clean-cc-stag2-s3p8 on
origin/clean-stagger2-latest = v17b 6bb6d0844; not yet in LATEST):
- **Old prod never ran k1_redist and face FCT together.** HG prod: face_fct on and time2_k1_redist=0 (the 0.1 on its
  command line is overridden by a later =0 on the same line). RSG prod: k1_redist=0.1 and no face FCT. With both on,
  k1_redist drains donors to exactly 0.1 E (c1 K1_E = -0.9 E_s) and FCT then cuts their inflow -> E_old < 0 that no face
  limiter fixes. Check your clean-cc keys for the same combination.
- **One-pass donor limiter** (`time2_face_fct_form = donor`, c298dfabd): a cell that stays above th V E^n even with its
  inflow cut keeps R = 1; every other cell takes R = (V(E+s) + I_safe - th V E^n)/P. Positive by construction, one R
  exchange, no Zalesak rounds, one flux per face (energy exact).
- **Cell sink as an outflow** (c9a918a83 + c1b9524b7): the cell's own exchange sink that its gas takes counts as one more
  outflow, scaled by the same R, bounded by the gas eint gain (dEN minus the KE of the momentum increment). Energy exact.
- Optional, as old prod did: scale the increment's face F by the same R (`time2_face_fct_f`, 4267f8d75; not run yet).
- HG MI300A 767342-346, v17b keys: base Zalesak 5.59 s/cycle, up to 775 rounds per call, 1.3e4 cells left Q<=0;
  donor + sink 4.81 s/cycle, 1 round per call, 0 CW, dt 97.5 held, E_old<0 per stage 163/22, 30/0, 12/0 after the
  restart cycle (the restart cycle itself has a K1 face-flux mismatch from U^n's F). r136cL: E_old<0 1906 -> 0, a=b=c bitwise.
- RSG A100 runs (pHG262 Raven) are next; small-deck gate + order test on CPU running.
Your 19o regime (E<0 in block-interior thin-per-cell cells from a large h div F) is the implicit side; the lagged donor +
damping + donor-implicit face source from our previous NOTE is what fixed it for us.
