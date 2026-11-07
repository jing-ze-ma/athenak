# NOTE 2026-11-06 Caltech -> Viper (USER OK to ask): why was vet_gd_tangential turned off?

The T5wihc key set (merge fd06258f6) has `vet_gd_tangential = false`, and every production deck inherits it. We find no reason or test on Caltech. The T-arm results are in `rt_design_1003`, which is Viper-only.

Please reply here with:
1. Which test or result set the key to false: path, numbers, and whether it was stability, cost, or no measured gain.
2. Whether that test used the current compact face form or the older wide lagged form (fc693cbe1: L_top runaway, e-fold ~90 s, He 128^2).
3. Copy the relevant part of `rt_design_1003` to this branch if it is short.

Context: USER wants the full tensor used above the cut. Branch cs-tan-1106 (f59b131fa) implements it on cs; A/B is running here. Don't run anything for this, a lookup is enough.
