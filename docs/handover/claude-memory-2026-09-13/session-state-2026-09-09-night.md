---
name: session-state-2026-09-09-night
description: START HERE next session (written 2026-09-09 ~22:45 CEST on orion) — what is running (prod11 chain, T15 prod12 gate), the one unfinished agent and its transcript, the uncommitted tree, and the decision tree for prod12
metadata:
  type: project
---

## Running (SLURM, orion; all under /orion/ptmp/jinma/Athenak/red_giant/)
- prod11 (lidded production): job 194624 running, 194625 queued (dependency). At 22:00 it
  was at t = 3.05e7, dt 30.1, healthy; bulk convection developed (rms v_r 0.5-1 km/s,
  ~12 cells around the equator at 0.95 R, photosphere at 1.097 R = right under the lid).
  Figures page: https://claude.ai/code/artifact/71df6ecc-d585-415c-b767-a9990c0c403f
  (scripts rglib.py/mkfigs.py were in the session scratchpad — gone; the agent's method is
  in the page's footer; docs/handover/scripts/dhjcs.py has the panel frames).
- T15_prod12gate: open top + sponge OFF (T12 config), tlim 3e6, launched ~22:35 by an
  Opus agent (job id in its transcript, see below, or `squeue -u jinma`). THE prod12 GATE:
  must reach 3e6 with no cold vertex cells. If it passes: prod12 = T6_nodust/rg.athinput
  + opac_tmin 2500 + sponge=false + the fixed binary (build a fresh pinned one from the
  current tree: it has the ghost-fill basis fix that 887c241e lacks), tlim 1e8, rst 5e6,
  bin 2e5, 6 nodes, chain like prod11/sub.sh.

## Unfinished agent (transcript survives; its LAST assistant message is the report)
/u/jinma/.claude/projects/-orion-u-jinma-ATHENAK-athenak/dbba3fb7-754f-4b17-8a81-74e65ff9f7ee/subagents/agent-a213074acdec93b96.jsonl
Tasks: (1) launch + monitor T15 to 1.5e6; (2) the HISTORY of the vertex chimney
(rho_vertex/rho_bulk at i=310 vs t through the T6..T11 open-top chain, in T12 sponge-off,
E3 wall, prod11). If the transcript has no final report, redo task 2 (the question is in
[[red-giant-vertex-chimney]]). Bin-dump trap: the first angular axis is REVERSED vs
mb_geometry; m=0's cube-vertex cell is (k,j)=(7,0), NOT (0,0).

## Decision tree for prod12
- T15 passes 3e6 AND the chimney does not deepen without the sponge -> launch prod12.
- T15 passes but the chimney deepens monotonically -> the open top has a vertex mass-loss
  defect; find it before a 1e8 run (the vertex cell is the only cell with two seam faces:
  flux_seam_cc.cpp 196-241; and the open ghost fill continues an already-evacuated column).
- T15 dies -> fall back to the lidded config (prod11 never dies) or damp only the 24
  vertex columns above the photosphere.

## Code state (repo branch polar-average-perf; origin at 40aef0ad)
PUSHED today: 5c0b98e4 WB restart cache, 32273ebb its note, 69919eae/79167f40/40aef0ad
docs/handover/NOTE-2026-09-09-cs-covariant-basis.md (viper's to-do list for the basis bugs).
UNCOMMITTED (the tree builds; build_rg_prod = 887c241e predates the ghost-fill fix):
red_giant.cpp (sponge CsKinetic fix, ghost-fill/IC/inner-pass basis fix, opac_tmin,
bg medium, vpert_rmin, MLT machinery incl. the chi-limiter removal to REVERT, rg_radbud,
open-ghost guard, sponge guard, rt_* switches), two_stream_rt.hpp (direct source + Newton,
ON by default, NOT A/B'd on the hot Jupiter -- hold), conduction rad_kappa_rmax/above,
hydro wb_rmax, driver.cpp every-cycle NaN checker, mesh/bvals cs_corner_poison debug flag
(default off), eos_dump.cpp argv. Commit in pieces when the user asks; the user asked
"did you push the fixes relevant to the hot Jupiter" -> answered: only the WB one, the RT
module needs a dhj A/B first.
Copied source tree with the ghost-fill fix + its build: /orion/ptmp/jinma/Athenak/src_vfix
(build_vfix); the main tree now has the same red_giant.cpp.

## Findings today, in order
[[red-giant-dust-opacity-kills-lidfree]] -> [[red-giant-restart-radial-ke-injection]] ->
[[red-giant-nan-check-hides-origin]] -> [[red-giant-sponge-cs-kinetic-energy-bug]] ->
[[cs-orthogonal-ke-audit-2026-09-09]] -> [[red-giant-vertex-chimney]] (the current state).

## ADDENDUM 2026-09-10 13:45 — where the prod12 thread is now
Read [[red-giant-explicit-conduction-explosion]] and the UPDATEs in [[red-giant-vertex-chimney]]
first. Gate = T19_blend30 (195144->195145) / T19_blend100 (195146->195147) to 3e6. An Opus
agent is polling them and, on a pass, LAUNCHES prod12 (fresh t=0, 1e8, 6 nodes x 3-job
chain, binary = src_vfix/build_vfix (ghost-fill fix, rad_dt_face/open_guard_fallback
default off), input = T19's + sponge=false + rad_tau 100/1000 (or 30/300), nan_check 10):
transcript .../884370bd-7441-4eee-9239-e5fc61707ea3/subagents/agent-a9a037df7e9660af4.jsonl.
If cut off: check `squeue -u jinma` for rgprod12; if absent and the T19s passed, do PART 2
of that agent's prompt by hand. Known blemishes accepted for prod12: 24 vertex chimney
columns (non-lethal, sponge off); explicit-branch stiffness moved below the front by the
blend. User has NOT been asked explicitly for the prod12 go -- the 09-09 plan said "if it
passes, prod12 = this config"; tell the user it was launched under that plan.

## ADDENDUM 2026-09-10 01:00 (session 293c860c) — the state after the hydro-trigger finding
Read [[red-giant-vertex-floor-hydro-trigger]], [[red-giant-molecular-opacity-knee-runaway]],
[[red-giant-prod11-died-3e7]] first. User asked (00:55) to FIX (1) the cube-vertex chart
asymmetry and (2) the energy-creating floor. Two Opus agents, each in its OWN source copy:
src_vtx/build_vtx -> bin_vtx (agent aa3e87412c077b925: onset dating, symmetric reproducer
with vpert off, bisection of flux_seam_cc two-seam-face cell / corner fill / resample end,
verify with the F0 config) and src_floor/build_floor -> bin_floor (agent a87c9fc892f825a2c:
hypothesis = the tabulated EOS clamps T at eos_logt_min=100 K and returns p(100 K) for a
cell whose e is far below -> pressure without energy -> Riemann creates energy; fix = p,
c_s scaled linearly in e below e_min, floors counted; verify on the F0 config with an A/B
switch). Running jobs: T22_tmin3200_r14x 195239(+195240) = the prod12 gate (clean at 3.26e6,
tlim 3.5e6), prod11_catch (agent abcfdfcf, reproducer of prod11's 3.2106e7 death, ~3 h).
Transcripts: .../293c860c-0e33-4e6c-bfe8-7387e9c93220/subagents/agent-<id>.jsonl.
