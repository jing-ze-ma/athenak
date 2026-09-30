---
name: noise-gate-equal-time
description: noise/A-B gates must stop every arm at the same PHYSICAL time (time/tlim), never by cycle count (nlim) -- dt differs ~5 % between arms and the steady deep drift then fakes a deviation
metadata:
  type: feedback
---
09-28 ckrev_0928: a HIP noise gate stopped by nlim showed ck-next "failing" 3-4.6x in the deep and a bizarre kick
dependence; arms had ended up to 320 s apart in physical time (dt set by the chaotic top, ~5 % spread) and the deep
drifts ~1e-6/s in T. With time/tlim the same gate PASSES (<= 1.05-1.44x vs 1e-16 members) and deep deviations drop 50x.
**Why:** discovered by the ck review worker; nearly led to a wrong revert recommendation.
**How to apply:** every noise / A-B / order gate brief: stop arms with time/tlim at an identical time; use a
round-off (1e-16) member ensemble as the yardstick; report arms' final times. Past nlim-based gates (e.g. Caltech's
ck-next 0.78x, ck_store_split 1.21x) may need a recheck if they matter.
