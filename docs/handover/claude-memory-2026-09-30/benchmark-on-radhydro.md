---
name: benchmark-on-radhydro
description: USER 09-26 - judge M1/VET defaults and speed levers on rad-hydro setups with MOVING gas (sph_wedge, He box), not on static radiation-only tests (T-S4 sph_atm)
metadata:
  type: feedback
---
User 09-26: "we want rad hydro so probably we should trust rad hydro test more."
**Why:** the static T-S4 wedge misled twice. Predictor off looked -6 % there but is +28 % on moving gas. mg_gc beat rbgs_fwd by 18 % there but is 8.7 % slower on the moving sph_wedge.
**How to apply:** time every lever and default on a moving-gas rad-hydro state: m1_test = sph_wedge (evolved restart, top sponge) or the He box. Static tests are for correctness only. Related: [[measure-cost-on-gpu]].
