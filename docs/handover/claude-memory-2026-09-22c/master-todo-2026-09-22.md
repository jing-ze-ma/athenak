---
name: master-todo-2026-09-22
description: MASTER PLAN agreed 2026-09-22 02:15 - order = wait for the running agents, CLEAN UP the code, then the three continuation lists (dhj RT, cs vertex/seam, cs dhj productions); M1 stays paused
metadata:
  type: project
---

User 09-22 ~02:15: "wait for all agents to finish, then clean up first, then proceed" (several lists).

## 0. In flight at writing (do NOT start new work before these land)
- CPU agent: IMPLICIT ck column solve (tests_ck_implicit/) + the repair of the ck_spherical=false path to
  the production plane-parallel form (regressed hybrid at HEAD, see [[dhj-improvement-list-0921]]).
- GPU jobs (all with ck_spherical = true unless noted): ck_sph_ab sphbeam 11931204 / prodbin 11931220
  (original binary, plane-parallel); recon plm_ctl/ppmx/wenoz 11931395-7; noang 11931398 (control =
  sph 11931203); wb-cache every1/every10 11931399/11931400 (apu, 10 h).
  READ: bench/*/README.md comparison plans. Results to record here when read.

## 1. CLEAN-UP (first, before any list) -> [[code-cleanup-pause-plan]]
Inventory every runtime switch (two-stream/ck, conduction, WB, rad_m1, dhj/box/RG pgens) with default,
production values, gate coverage; user decides per switch: ABANDON / DEFAULT ON / KEEP. Small commits,
regenerate bitwise references when a default changes. Also: repair or remove the HEAD hybrid; stale
bench NOTES.md; restart-format versions; docs/handover update; push to the fork after user OK.

## 2. LIST A - dhj radiative transfer (after clean-up)  -> [[dhj-improvement-list-0921]]
a. Read the GPU A/B: production (prodbin) vs sph vs sphbeam: dt (sph 12.5 vs prod 20.8 s: WHY), floors,
   T(p) day/night at 1e-6..1e-3 bar, cost. Then decide defaults (ck_spherical, ck_beam_sph on?).
b. Implicit ck column solve: gates, cost; fold conduction (removes the 10 bar split); then default?
c. rad_angular = false for dhj if the noang arm is inert (committed switch 213095e4).
d. Beam power reaching the ck cut is dropped (0.003 %): fine, documented.
e. Long-term: short characteristics per band/g-point (heating, beam, VET tensor). M1 stays PAUSED
   ([[m1-paused-back-to-dhj]]).

## 3. LIST B - cs vertex / seam accuracy  (tests_cs_regions/ is the gate)
a. Angular ppmx or wenoz instead of plm (vertex error x0.31 L1, x0.15 Linf, cost 1.28x); ppm4 is bad.
   GPU recon arms decide for the real atmosphere.
b. Make the RADIAL (x1) reconstruction selectable on cs (always GridPLM today).
c. REAL cube-vertex fill (sampling the third panel; prototype beat extrapolation 3.4-5.5x) - the one
   real build; plus the corner EMF at the 3-valent vertex (field-loop vertex L1 unaffected by any
   reconstruction -> CT/EMF suspect).
d. cs_wellbalanced_src: does NOT help the vertex (8 % worse), unreadable from <hydro>: fix or abandon.
e. Receiver-side seam resample (3rd-order block-end residual): low priority.
f. Closed: seams (2nd order), fc halo (no defect), cc corner slots (fixed ec0a9724), angular momentum
   (cs not worse than sp, tests_cs_angmom/).

## 4. LIST C - cs dhj production practice
a. wb_cache_every = 1 if WB is used (GPU arms test the rot-12.3 death); cs_mhd_prod3_wb ran WITHOUT
   cs_wellbalanced_src (its NOTES.md is wrong).
b. dt_min 1e-2 -> 1e-3.
c. <mhd>/fofc was never on; dhj.log fofc column is meaningless; use eos_dfloor/efloor/tfloor/eos_fail.
d. Restarts: MHD + cs caches now bitwise (0eb7e2c4, 544cb36b); pgen state (dump clocks, WB cache
   flags) not audited; cells flooring in the restart cycle.
e. Floors at the night-side vacuum top (dfloor vs WB background, keep-temperature).
f. Next production: new input with the chosen defaults; chain from scratch or from rot 283?
   (user decision). bench/cs_mhd_prod3/rst = 46 GB (keep a subset? user).
g. Physics questions parked: shallow jet / deep westward flow (sp agrees), day-night contrast at low p.

## 5. Also parked
- He4 global star (he4-session-state-2026-09-21), B-star/He boxes: see their notes.
- Upstream merge (72 commits behind origin/main; 6 overlapping files).
- origin remote URL embeds a token (user action).

## Update 09-22 02:40: implicit ck agent DONE (committed after 76359975). The HEAD hybrid is REPAIRED
(ck_spherical=false = production plane-parallel divergence; sweep-to-gas gap -8/-4 %, not -36 %).
ck_implicit: gap 1e-9..3e-8 in all combos, but 5.1-5.4x RT cost and DIVERGES at 10-100x dt (thin-cell
linearisation; needs the grey mode-3 two-level split rt_impl_tau_min) -> NOT production-ready; list A item b.
ALL CPU AGENTS DONE. Remaining in flight: the GPU jobs only. Next = CLEAN-UP (section 1).
