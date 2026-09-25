---
name: cs-mhd-prod4-run
description: cs_mhd_prod4 LAUNCHED 2026-09-22 ~20:45 - the cubed-sphere dhj production with the corrected RT (spherical tm sweep + pseudo-spherical beam), polytropic WB, from scratch; jobs 11941995-8; what to check at rot 0.5; binary swap at a link boundary APPROVED once the speedups pass (see body)
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

**USER APPROVED 09-22 ~21:45: BINARY SWAP AT A LINK BOUNDARY.** When the IMP/LEG template-tag patch and the
rank-imbalance remedy (leading/trailing split) each pass their own gate (bitwise dumps vs HEAD + GPU
speedup), build ONE clean-snapshot binary with the passing change(s) (HEAD + patches, same recipe as
prod4: HIP GFX942 MPI Release, PROBLEM=deep_hot_jupiter_rt; if the split is a runtime switch, turning it
on in the prod4 input is part of the swap). Before swapping: restart a prod4 restart file on apudev with
old and new binary for a few hundred cycles -> dumps/restarts must be bitwise identical. Then replace
bench/cs_mhd_prod4/athena BETWEEN links (never while a link runs; keep the old one as athena.395db5bc),
record md5s + reason in its NOTES.md. Only one of the two passing -> swap with that one alone.

**09-22 ~22:20 USER: prod4 HELD (scontrol hold 11941995) until the rank-imbalance agent reports. Release with scontrol release 11941995 after the decision (swap-in the split binary first if it passes).**
**09-22 ~21:45 SWAP GATE (bench/prod4_swap_0922, job 11942248): payloads BITWISE in all arms (old/new, 12/12 and
10/14; header differs only by the lb default line), but from a COLD START 10/14 is SLOWER: cpu 24.73/24.88 s vs
22.07/22.01 (old 12/12), new 12/12 22.25. The rot-283 measurement (tests_rank_imbalance) said +14.5 %: the best
split is state-dependent. prod4 RELEASED with the OLD binary. Swap binary kept: athena.swap md5 e5cfbc99.
Next: re-measure 12/12 vs 10/14 from prod4's OWN restart once it has developed (e.g. rot 10+); swap at a link
boundary only if faster there.**

**09-23 00:45 cs_wellbalanced_src set FALSE in the run-dir input before the first link started (user: drop it if it does not help).** Min beta at prod3 rot 283 = 1.22 (>> ~2e-3 crossover); hurts hydro accuracy. Both runs changed together (twins). Smoke of prod4 input with it off: job 11943851 (bench/prod4_wboff_smoke) -- check its result before 02:36.

**09-23 ~03:30 SPEED-UP AVAILABLE, needs the user: `export HSA_NO_SCRATCH_RECLAIM=1` in submit.sh** cuts prod4 GPU time 22.5 -> 17.0 s per 300 cycles (1.32x) with bitwise identical output (tests_m1/runs_3k_gpu3d/README_HALO.md, job 11944726). Links 2-4 are already queued and sbatch snapshotted their scripts, so applying it = edit submit.sh + cancel and resubmit the pending links (dependency on the running link). Same for cs_hyd4_prod.

**09-23 ~03:40 APPLIED (user):** HSA_NO_SCRATCH_RECLAIM=1 in submit.sh; pending links resubmitted (prod4 11944873-5, hyd4 11944876-8); running first links (11941995, 11943719) unchanged.

**09-23 03:25 ROT ~1.2 CHECK (once):** prod4 tracks prod3's cold start: ME 8.19e32 at rot 1.1 vs prod3 8.04e32 at rot 1.2 (winding-up growth, not the -5 %/rot decay of the developed state); dt 3.8 s at rot 1.26 vs prod3 4.35 s at 1.2 (prod3 then rose to 10.4 s over rot 1-5, 16.6 s over 5-20); eos_fail 0; per-interval dfloor 4.0e8, tfloor 3.2e7, efloor 1.6e7 -> 9.6e6 (falling, plateauing as expected). Pace now ~0.66 rot/h (low early dt); 4 links will not reach 283 rot at this pace -- re-assess after rot ~20 when dt should reach ~16 s.

**STOPPED 09-25 03:30 (user: cancel everything in prod).** All links cancelled; newest restarts: prod4 dhj.00179.rst (~rot 89.5), hyd4 dhj.00217.rst (~rot 108.5). Superseded by the new physics (ck_nquad 2, ck to the wall, ck RCE IC). Can be continued from these restarts with the old binary if ever needed.
