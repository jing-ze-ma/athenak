---
name: session-state-2026-09-09-evening
description: START HERE next session (written 2026-09-09 ~18:10 CEST on orion) — what is running, what the two unfinished Opus agents were doing and where their transcripts are, what is uncommitted, and the decision tree for prod12
metadata:
  type: project
---

## What is running (SLURM, orion)
- prod11 (lidded production, /orion/ptmp/jinma/Athenak/red_giant/prod11): job 194515 at
  t = 8.56e6, dt 30.3, healthy; convective KE_r grew x2e4 since 3e6 (bulk convection igniting),
  L_out/L plateau 0.63-0.68, L_cut/L 0.20-0.24. Chain 194624 -> 194625 restarts with
  `hydro/wb_cache_every=1` (protects against [[red-giant-restart-radial-ke-injection]]).
  Its vertex cells are 5-8% cold (sponge bug, sub-lethal there).
- T11_spongefix (job 194663): FAILED at 8.72e5 -- same vertex collapse as T7 (see the UPDATE in
  [[red-giant-sponge-cs-kinetic-energy-bug]]). prod12 is NOT unblocked. Original plan text: Restart from t = 1.0e6 (T9's rst) of the grains-off, lid-free, open-top star with
  the sponge fix (binary 887c241e), nan_check every cycle, tlim 3e6, bin every 2e4.
  Must: cross 1.04e6 (old death), keep vertex cells healthy (e/rho ~1.8e12, v_h <= bulk;
  sick = e/rho 1.42e7, 40 km/s), pass 1.5e6, reach 3e6. Check with the vertex script idea in
  [[red-giant-sponge-cs-kinetic-energy-bug]] (blocks with mb_logical lx2,lx3 in {0,3}, local
  cell k=0|7, j=0|7). If it passes: prod12 = this config (T6_nodust/rg.athinput + opac_tmin
  2500 + fixed binary), tlim 1e8, rst 5e6, bin 2e5, 6 nodes, chain like prod11/sub.sh.
- T10_spgguard (194653, agent's control: sponge positivity guard only, binary ac2b319c,
  tlim 1.06e6): tells whether the guard alone masks the drain. Read its out.txt.

## NEXT STEPS (in order)
1. Poison test of FillPanelCornersCC corner ghosts (1-cycle job): does any ACTIVE cell read them?
2. Sponge-off rerun from T11's 8e5 rst (T11_spongefix/rst/) to fully exonerate the sponge.
3. Read the vertex-adjacent edge-ghost resample at the seam end (bvals_cc seam resample / flux_seam_cc).
4. Apply the audit's ghost-fill fixes anyway ([[cs-orthogonal-ke-audit-2026-09-09]]).
5. Fallback if the vertex cannot be fixed quickly: keep the top pressure-supported (lidded
   prod11 never collapses) or damp only the 24 vertex columns above the photosphere.

## Finished agents (their transcripts survive on disk; the LAST assistant message is the report)
- NaN-origin agent (DONE; its A/B result is summarized above): /u/jinma/.claude/projects/-orion-u-jinma-ATHENAK-athenak/375b2275-3f04-4aca-aa3a-f35153b1627f/subagents/agent-a9abaadd14ebee186.jsonl
  (was asked to run/monitor T11 and report vertex health + the task-2 extremes table; it
  also improved the NaN checker in src/driver/driver.cpp -- active vs ghost counts, neighbour
  state, block->panel map -- UNCOMMITTED).
- cs orthogonal-KE audit: DONE, summarized in [[cs-orthogonal-ke-audit-2026-09-09]] (apply its fixes before prod12).
  (asked for a per-site table: prolong_prims.cpp, ideal_c2p_*.hpp, pgen.cpp, red_giant.cpp's
  7 other sites, deep_hot_jupiter.cpp sponge/drag, viscosity/conduction, newdt, fofc). If the
  transcript has no final table, rerun the audit (prompt is in this session's transcript).

## Code state (repo /orion/u/jinma/ATHENAK/athenak, branch polar-average-perf)
- COMMITTED+PUSHED today: 5c0b98e4 (WB cache rebuilt on first call after restart -- affects
  every WB restart incl. dhj; viper must pull), 32273ebb (docs/handover/NOTE-2026-09-09-wb-restart-bug.md).
- UNCOMMITTED (build_rg_prod = 887c241e has all of it): red_giant.cpp (bg medium, opac_floor
  table extension, opac_tmin, vpert_rmin, MLT machinery incl. the chi-limiter removal that
  should be REVERTED, rg_radbud, open-ghost guard + column fallback [make it a clamped
  continuation before production], sponge positivity guard, CsKinetic + sponge fix, rt_*
  switches), two_stream_rt.hpp (direct source, closed-form+Newton, rt_explicit, top slot from
  cell ie, Newton positivity guard), hydro wb_rmax, conduction rad_kappa_rmax/above,
  driver.cpp NaN checker, eos_dump.cpp argv range. Commit in pieces; sponge fix + opac_tmin +
  NaN checker are the clean ones.

## Findings today, in order (each has its own note)
[[red-giant-dust-opacity-kills-lidfree]] (deaths at 2-4e5) -> [[red-giant-restart-radial-ke-injection]]
(x730 kick, fixed) -> [[red-giant-nan-check-hides-origin]] (why no precursor) ->
[[red-giant-sponge-cs-kinetic-energy-bug]] (deaths at ~1e6, fixed, verifying). Results page:
https://claude.ai/code/artifact/2186f16c-766f-470d-86c7-dec2084bd1c3 (sections 1-18).
User preference: I diagnose and verify, Opus agents do the legwork; but the user also said
"you probably need to find the bug yourself instead of handing it to opus" when rounds
multiplied -- do the thinking first, delegate the mechanical part.
