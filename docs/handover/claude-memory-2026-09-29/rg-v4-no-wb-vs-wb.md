---
name: rg-v4-no-wb-vs-wb
description: 2026-09-13 RG_v4 (red giant from scratch, NO WB, exact-flux wall d2788159, FOFC, EOS clamp) vs RG_fofc (WB) over 0-9e5 s. Conservation is better without WB (tot-E drift 0 in 6 digits, mass exact, no dt collapse, floors/FOFC below the WB run) BUT the velocity field is worse: a STEADY spurious mean radial flow from the first dump (i=0 -3.1e4 cm/s, i=1-4 +6e3 to +1e4, odd-even, frozen to 0.7% over 8e5 s; +590 cm/s outward mean at i=20 = 3.6x the rms) vs 950 cm/s at i=0 and ~0 in the convective zone with WB = the hydrostatic truncation residual WB cancels. Base convection still switches on (i=10 e-fold 3.0e5 vs 2.5e5 s, rms 567 vs 480 at 9e5) but is confined to i=4-8 with a sharp peak at i=7 (724 cm/s), marginally STABLE at i=9.5, and COUNTERGRADIENT at i=12-20 (C(v_r,T') -0.45, C(v_r,rho') +0.85) where the WB run is buoyant (+0.9). Neither run crosses 1e4 cm/s by 9e5 (WB run: i=10 at 2.0e6, i=5 at 2.5e6)
metadata:
  type: project
---
Scripts bench/RG_v4/analysis/cmp_{vr,eos,hst,cross}.py -> *.txt (mirrors of RG_fofc_long3/analysis/partb.py,
parta_realeos.py). Caveat: RG_v4 also carries the newer binary (EOS clamp, corner limiter, exact-flux wall), so the
floor counters and part of the base difference are confounded with WB per se. Decision for the user: WB (small
cancelled residual, 1e-6 E offset, wb_rmin/wb_rmax bounds) vs no WB (exact conservation, 30x larger steady wall/envelope
mean flow that forces the base). Middle ground not yet run: WB with the exact-flux wall (RG_fofc_long3 predates
d2788159). See [[red-giant-inner-mode-is-convection]], [[red-giant-3e6-long2-vs-long3]], [[rg-box-sweep]].

FOFC burst 1.46-1.54e6 s (rg.log fofc 0 -> 3.2e5/interval, dfloor 2.3e6): analysis/fofc_where.py -> 100 % of flags and floored
cells at i>=474 of 480 (open top, above the photosphere i~338), zero at the base/envelope, anti-correlated with panel edges;
same character as RG_fofc_long3's burst at 1.06-1.14e6 (peak 5.8e6, then decayed), 20-50x weaker and 4e5 s later. Corona-only,
the open-top/corona-drain physics decision, not a numerics problem.
