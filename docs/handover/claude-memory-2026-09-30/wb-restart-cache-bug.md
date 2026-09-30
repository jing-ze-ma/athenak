---
name: wb-restart-cache-bug
description: Orion 09-09 - the well-balanced background cache was all zeros for up to wb_cache_every-1 cycles after EVERY restart (FIXED 5c0b98e4); a coherent radial kick per chain restart in every pre-fix WB run
metadata:
  type: project
---

`wbq0` was rebuilt only when `ncycle % wb_cache_every == 0`, and ncycle is restored from the
restart while the cache is freshly zeroed, so 9 of 10 restarts ran up to 9 cycles with NO
background: a free-fall kick g*(wb_cache_every-1)*dt, spherically coherent. Fix 5c0b98e4
(`wb_cache_built` flag, zero cost). Stopgap for a binary that must not change:
`hydro/wb_cache_every=1` (and mhd/) on the restart line. Signature in .hst: 1-KE jump in
the first rows after a restart. Details: docs/handover/NOTE-2026-09-09-wb-restart-bug.md.
Orion also showed the CPU restart path is bit-identical across rank layouts, so viper's
5e-6 restart divergence is GPU-specific. NOT the cause of [[cs-vertex-dt-collapse-0907-defaults]].
