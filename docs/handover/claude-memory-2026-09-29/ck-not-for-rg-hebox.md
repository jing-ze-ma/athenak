---
name: ck-not-for-rg-hebox
description: user 09-27: red giant and He box do not need correlated-k; no ck tests/defaults for those pgens
metadata:
  type: feedback
---
Correlated-k (rt_ck, ck_nquad, ck_implicit ...) is for the hot Jupiter (dhj) only. Red giant and He box
(box_convection) do not need it.
**Why:** user 09-27 night, when a ck_nquad 2 default test for red_giant/box_convection was proposed.
**How to apply:** do not propose ck switches, tests or defaults for red_giant / box_convection; their radiation is
the grey two-stream / M1 path.
