---
name: red-giant-rt-bface-switch-and-head-regression
description: "09-11 05:50: problem/rt_bface (neighbour-Planck BFace fix) defaults FALSE in pin8/pin9 (pin7 and earlier had it UNCONDITIONAL) -> every run input must set rt_bface=true; HEAD regression: solar_convection + picket-fence dhj bit-identical, correlated-k dhj differs ONLY through the 4 new-by-default RT switches (rt_src_direct/rt_newton/rt_rescue_eq=true, rt_semi_lin=false) + the unconditional band top-slot clamp; user decides whether to flip them"
metadata:
  type: project
---
Verified: `strings athena_pin7 | grep -c rt_bface` = 0, pin8/pin9 = 2. The R2/I6 inputs never
mention rt_bface, so any pin8+ run without `problem/rt_bface = true` runs the OLD Planck
endpoint (the cold-collapse bug of [[red-giant-rt-neighbour-planck-bug]]). R3_scan job 196156
(pin9, first launch 05:35) had to be relaunched with the flag; Agent A was told. The repo input
inputs/hydro/red_giant_cs.athinput now sets it.
Regression harness: /orion/ptmp/jinma/Athenak/regress/ (b_{wt,hd,mod*}_{sc,dhj} builds,
run_* dirs, cmp.py); worktree removed. Full inventory of gated / un-gated / diagnostic hunks:
docs/handover/NOTE-2026-09-11-red-giant-switches.md (untracked).
Recommendation (not yet decided by the user): flip the 4 RT defaults to old behaviour, gate
the top-slot clamp, set all five in the red-giant input, re-check bit-identity both ways.

## DONE 09-11 06:30: RT defaults flipped to OLD behaviour (verified diff + bit-identity)
two_stream_rt.hpp: rt_src_direct/rt_newton/rt_rescue_eq default false, rt_semi_lin true, NEW
problem/rt_top_clamp (default false) gates the band top-slot clamp at 4 sites; red_giant.cpp both
configure blocks updated. Gates: ck dhj vs HEAD bit-identical (regress/run_ck_new vs run_ck_hd);
red giant new build + opt-ins vs pin9 bit-identical (run_rg_new vs run_rg_pin9, 100 cycles).
CONSEQUENCE: every red-giant run input from pin10 on MUST carry the SIX opt-ins (rt_bface,
rt_src_direct, rt_newton, rt_rescue_eq = true, rt_semi_lin = false, rt_top_clamp = true) --
inputs/hydro/red_giant_cs.athinput has them; the ptmp run inputs (I6_fixes/rg.athinput etc.) do
NOT. pin9 hard-wires the new solver and is unaffected. Docs note §3 verdict updated.
