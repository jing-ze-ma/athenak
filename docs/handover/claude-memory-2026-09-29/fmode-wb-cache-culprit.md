---
name: fmode-wb-cache-culprit
description: "f-MODE GROWTH FOUND 09-16 17:20 = the STALE well-balanced background cache (wb_cache_every = 10 in every production input); fresh cache (0/1) or WB off kills it (gamma -0.07/-0.04/+0.05 vs +0.86, 6-7 sigma); dt scan on the merged binary reproduced the old line first"
metadata: 
  node_type: memory
  type: project
  originSessionId: 906729a4-d9ad-4c83-98b9-9ada51de1fc9
  modified: 2026-09-16T15:13:27.405Z
---

Small B-star box 336x64x64, cfl 0.3, merged binary b775a270 (md5 887d1471), bench/fmode_0916b/wb/ (jobs 11728412-4, RESULTS.txt, fit_wb.png):
  N03 reference wb_cache_every=10: gamma +0.864 +- 0.108 per turnover (reproduces old F03 exactly)
  W0 WB off (wb_x1=false + wellbalance_dynamic=false): -0.067 +- 0.090
  WD wb_cache_every=0 (rebuild every stage): -0.039 +- 0.081
  WC wb_cache_every=1 (rebuild every cycle): +0.054 +- 0.082
  Cost of a fresh cache: not measurable (81 vs 86 ms/cycle, node noise).
Reading: the hydrostatic WB reference (hydro_fluxes.cpp ~170, hydro.hpp wb_cache_every) is rebuilt only every 10 cycles, so in a flowing background gravity is balanced against a background up to 10 dt (~15 s) old: a lagged O(dt) forward-Euler-like kick on the surface gravity wave = gamma = 1.31 omega^2 dt. Explains: needs convection (static column exact), even under RT operator exchange, survives Strang/ImEx of the radiation. Note wb_x1=false alone FATALs ("wb_direction not set"); wellbalance_dynamic=false alone == WB off (every site gates on both).
dt scan rerun same day (bench/fmode_0916b, N03/N015/N0075): new binary reproduces the old line (0.86, 0.42; intercept -0.02 +- 0.19) BEFORE the cache fix.
Production inputs all carry wb_cache_every = 10: bstar prod_w4_p, He box_w4_in / _in2. Fix = wb_cache_every = 1 (or 0).
CONFIRMED 09-16 19:00 (bench/fmode_0916b/wb/fit_wc.py, arms WC/WC015/WC0075, cache 1): gamma +0.054+-0.082 / +0.004+-0.068 / -0.033+-0.057 at dt 1.50/0.75/0.37 s; slope 0.08 per s (was 0.61), intercept -0.06+-0.08 = the physical damping. The dt dependence is GONE. Productions prod_w5/box_w5 relaunched from t=0 with wb_cache_every = 0 (see MEMORY.md CURRENT STATE). Still open: the production-box mode amplitude must be seen to decay in rt_surface.bin; port wb_cache_every = 0 to rg/dhj. See [[fmode-dt-taper-test]], [[rt-source-dt-forcing]], [[cs-wb-source-cached]].
