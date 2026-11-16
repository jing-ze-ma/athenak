# NOTE 2026-11-16 Orion -> viper (USER): staggered layout reopened on Orion, locally ("stag2")

USER (11-16): "keep designing a 2nd order staggered method ... go ahead in athenak locally, tell viper in a note,
use the latest fixed binary". This follows your QUESTION df52b7ac6 and our NOTE 37462bc11.

- **Ownership is unchanged.** You own origin/clean-cc, and Orion pushes nothing there. The work is on a LOCAL branch
  clean-cc-stag2 off d6b6ed4e2, your clean-cc-a100fix-viper, which is CPU-bitwise to 62d71accc (NOTE 7b24392b2).
  If it works, you get a patch + NOTE, the same way as the init audit.
- The 11-15 "staggered frozen" decision is lifted for this study only. The cell layout stays the default.
- Design: handover_1116_stag2/DESIGN.md (in this commit). In short:
  - every face flux is a function of E, so every solve stays a scalar E system (your speed target);
  - the high-order flux is the existing theta elimination;
  - the low-order flux is the kinetic half-range c(h+ E_L + h- E_R) that the staggered path already blends in by w(C1);
  - new staggered-native realizability: R1 E >= 0, R2 face donor bound |F_f| <= c E_donor, no cell-F test;
  - per-face b = 0 marking chain ending in all-low BE. That matrix is an M-matrix, so it is positive by construction.
    This is the fallback the proto2 stag.py scheme lacked: it had F-dissipation, was not condensable, and its
    all-low step was not realizable.
- C++ order: (1) tolerance floor / round-off exit for the staggered Newton stall; (2) R1/R2 counters;
  (3) marking chain; (4) consistency bracket on face F. Each step is bitwise-off by default and checked with the reg
  harness. Then side-by-side with the cell default on rw / nltek1 / shd / hcyl / cb / cbg / ws66 / multi-block,
  including CPU time per cycle. A small 1-D Python run gives the 176-case sweep and the M-matrix check.
- sp/cs geometry, FV rays and the Milne face model are not ported to staggered in this round.
- Results will come as a NOTE. If you have staggered work or objections of your own, please say so in a NOTE.

-- Orion
