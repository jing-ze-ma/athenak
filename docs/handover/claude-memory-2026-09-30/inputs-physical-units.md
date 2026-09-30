---
name: inputs-physical-units
description: User 09-29: input parameters must be in physical units with explicit conversion; bbot was Heaviside-Lorentz code units (3 G = 0.846) - "3 G" in old notes (prod4) was really 10.6 G
metadata:
  type: feedback
---
User 09-29: "why is bbot the code unit? isn't that confusing? ... make sure all the input parameters are in physical unit
and properly converted". Found: dhj problem/bbot is Heaviside-Lorentz code units (B_G = bbot * sqrt(4 pi)); prod4
"bbot = 3" = 10.6 G, not 3 G. Audit + physical-unit keys (bbot_gauss etc.) delegated 09-29 ~22:15.
**Why:** code-unit inputs silently mislabel physics (a 3.5x field error went unnoticed for days).
**How to apply:** every new input key takes a physical unit named in the key or comment and the code prints
"value [code] = value [unit]" at startup; when quoting a field or any dimensional input, state the unit system; treat
old "bbot 3 G" statements as code units (10.6 G).

**MERGED 09-29 ~23:00:** units-physical-inputs on rt-integration (9e0b204a, linear; tree = gated merge ede1950f) + tfloor_kelvin = 5000 in He inputs (354a88b7), pushed. Gates: dhj hydro, MHD old key, bbot_gauss = bbot bitwise, He box bitwise (floor never binds), 9 rad/rad_m1 CPU tests pass. Use problem/bbot_gauss from now on.
