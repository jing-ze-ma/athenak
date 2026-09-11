---
name: session-state-2026-09-10-early
description: START HERE next session (written 2026-09-10 01:25 CEST on orion, session 293c860c) — what is running (FL1/FL2/FL3 floor A/B, V4 chart-free seed, prod11_catch reproducer), what to read first, the source copies and pinned binaries, the prod12 recipe, and the open questions
metadata:
  type: project
---

## Read first, in this order
[[red-giant-vertex-floor-hydro-trigger]] (the hydro trigger caught in the act) ->
[[red-giant-molecular-opacity-knee-runaway]] (what cools the cells; opac_tmin=3200 gate
PASSED 2.0e6->3.5e6 on the restart branch, fails from scratch) ->
[[red-giant-floor-energy-creation-fix]] (the fix; FL1 passes the death) ->
[[red-giant-vpert-seed-chart-imprint]] (vertex asymmetry = the seed, not the seams) ->
[[red-giant-prod11-died-3e7]] (the lidded production is dead at 3.2106e7, reproducer running).
The 09-09 explicit-conduction/rad_dt_face story is REFUTED; keep rad_dt_face as a diagnostic.

## Jobs running at 01:25 (all under /orion/ptmp/jinma/Athenak/red_giant/<dir>/out.txt)
- FL1_floorfix 195284 (both floor switches, 2 nodes, tlim 1.56e6): at 1.471e6 clean.
- FL2_ekin_only 195292 and FL3_consistent_only 195293 (one switch each, same config): were
  at 1.41e6; the baseline dies at 1.43114e6 -> check `grep -a COLLAPSE out.txt` and the
  efloor_de column in the event log (rg.*.log) -> which switch carries the result.
- V4_prod12cart 195275 (2 nodes, from scratch, open top, sponge off, vpert_cart=true, tlim
  3e6, ~4 h wall): does the vertex chimney / any death appear? Baselines: T6 open chain
  noise ratio > 1.5 at 2e5, first collapsed cells 4.4e5, deaths 7.7e5-1.04e6. Use
  _analysis_0910/task_floor.py (census2_run, vertex_copies) and task_noise.py on its dumps.
- prod11_catch 195261 (6 nodes, nan_check 1, rst every 2.5e5 in rst/, bin every 2e4 from
  3.19e7): the FATAL block expected ~03:35; classify seeds with
  _analysis_0910/task_G4_seed.py; then bisect from the nearest rst (sponge=0, opac_tmin).
The agents' transcripts (their last message is the report): .../projects/-orion-u-jinma-
ATHENAK-athenak/293c860c-0e33-4e6c-bfe8-7387e9c93220/subagents/agent-{a87c9fc892f825a2c
floor, aa3e87412c077b925 seed, abcfdfcf04aac973f prod11, ae9c52e5b4c0c5944 catch}.jsonl.

## Source copies and binaries (main repo UNCHANGED tonight; nothing committed)
- src_vfix/build_vfix: rad_dt_face + open_guard_fallback + rt_dump_rank (diagnostics).
- src_floor/build_floor -> bin_floor/athena + bin_floor/floorfix.diff (the floor fix: eos.cpp/
  .hpp, eos_table.*, general_c2p_hyd.hpp, general_hyd.cpp, coordinates.cpp, prolong_prims.cpp,
  mesh.hpp, eventlog.cpp; switches default OFF). Build with -D CMAKE_CXX_COMPILER=$(which
  mpicxx) after module load gcc/13 openmpi/4.1.
- src_vtx/build_vtx -> bin_vtx/athena (problem/vpert_cart, default off).
TO DO when merging: port floorfix.diff + vpert_cart into the main tree (they touch disjoint
files), add ThermoAt's sub-floor branch, run the style check, commit in pieces.

## prod12 recipe (not launched; user decides)
Binary = main tree + floor fix (both switches ON) + vpert_cart=true; open top, sponge off,
grains off (opac_tmin 2500) to ~1.4e6, then opac_tmin 3200 (a restart with the clamp), OR
test whether V4's chart-free seed + floor fix survives from scratch with 3200 directly;
nan_check_cycles <= 10, rst every 2.5e5, no 42 GB logs (UCX_LOG_LEVEL=error in sub.sh).

## Open questions
1. The momentum kick at the floored vertex cell (2.6e4x any force; agent a87c9fc was asked
   to decompose it from D1's fine dumps) -- WB source? Riemann with c_s 2100 vs 5e5?
2. Why open-top atmosphere cells see so little J (net cooling 15-100% of emission at
   tau ~0.3 above an 8000 K front) -- 1-D columns across a corrugated front?
3. prod11's 3.2106e7 death (reproducer).
4. The early open-top cooling to 2000-2500 K in the first 1e5 s (all open runs).
