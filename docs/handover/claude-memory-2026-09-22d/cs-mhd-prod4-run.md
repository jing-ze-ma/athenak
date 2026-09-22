---
name: cs-mhd-prod4-run
description: cs_mhd_prod4 LAUNCHED 2026-09-22 ~20:45 - the cubed-sphere dhj production with the corrected RT (spherical tm sweep + pseudo-spherical beam), polytropic WB, from scratch; jobs 11941995-8; what to check at rot 0.5
metadata:
  type: project
---

bench/cs_mhd_prod4, launched 2026-09-22 ~20:45 by the user's order. Binary md5 f95130b2..., clean snapshot
of 395db5bc (HIP GFX942 MPI Release). Jobs 11941995 -> 996 -> 997 -> 998 (apu, 2 GPUs / 1 node, 24 h
links, prod3's chain/STOP logic). FROM SCRATCH, tlim 8.64e7 (283 rot), outputs as prod3.
= prod3 + exactly the nine documented changes: ck_spherical + ck_beam_sph + ck_sweep_form = 1 (tm),
polytropic WB (wellbalance_dynamic, wb_x1, wb_cache_every 0, cs_wellbalanced_src), rad_angular false,
dt_min 1e-3, max_eta 5e12, pfloor 1e-5 barye; dead keys cs_gs07_emf / cs_full_rotation removed.
Semi-implicit gas coupling (cold-start split-lag gate passed, ratio 1.06). NOTES.md in the run dir.
**How to apply / at the first restart (rot 0.5) and on every look:** dhj.mhd.hst: dt should fall from the
cold-start CFL to ~13-20 s and stay CFL-limited (max_eta 5e12 keeps the Ohmic cap off); re-fit
d ln(1-KE)/dt for a slow dt-dependent mode; ME decay should be ~-5 %/rot (prod3), NOT -35 %; dhj.log:
fofc column is meaningless (off), watch eos_efloor (nonzero from cold start, 8.8e6/0.17 rot in the gate,
harmless but should plateau), tfloor/dfloor fractions ~0.5 % / 2 %; night-side T(p) at 1e-6..1e-3 bar
vs prod3 (expect the spherical+beam values, within ~50 K of prod3 at the top). Deaths at rot 11-12 in
old WB arms were NOT the cache; if this one dies there, the MHD dt-cell diagnostic names the cell.
Do not chain further or relaunch without the user. See [[dhj-improvement-list-0921]],
[[master-todo-2026-09-22]].
