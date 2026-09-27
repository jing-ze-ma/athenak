---
name: index-old-session-states
description: Older per-session state notes (orion) and the superseded trailing state paragraphs from MEMORY.md.
metadata:
  type: reference
---

## Session states (orion)

- [Session state 2026-08-16](session-state-2026-08-16.md) — superseded: four unpushed commits on general-eos; the open thread was the solar atmospheric runaway
- [Session state 2026-08-17](session-state-2026-08-17.md) — previous entry point (origin/general-eos = d1c289dd)
- [Session state 2026-08-18](session-state-2026-08-18.md) — superseded entry point: no code changed today; all cost questions answered; the ONE open thread is the viper dt discrepancy
- [Session state 2026-08-24](session-state-2026-08-24.md) — superseded: the ckrepro reproducers ran clean; hst is unusable here (cadence + Cartesian volume)
- [Session state 2026-08-25](session-state-2026-08-25.md) — superseded: job 190939 was left running that day; the viper doc's blow-up does NOT reproduce on orion CPU; input pushed as cd2e8815
- [Session state 2026-09-04](session-state-2026-09-04.md) — START HERE: artifact PUBLISHED with the time player + ray-marched convection plumes; nothing pending
- [Session state 2026-09-07 (orion)](session-state-2026-09-07-orion.md) — START HERE: pulled viper's handover (5e9b1273); read docs/handover/HANDOVER-2026-09-07.md; steps 1-2, the determinism test and the EOS table all DONE; left: GPU build with d3d74f2b (needs viper)
- [Session state 2026-09-09 evening](session-state-2026-09-09-evening.md) — superseded by the NIGHT note; keeps the T11/T10 details and the earlier agent transcripts
- [Session state 2026-09-09 night](session-state-2026-09-09-night.md) — superseded by the 09-10 EARLY note: prod12 gate FAILED (see the explosion note), nothing running except prod11 chain 194625; the user must choose: bigger domain, stay lidded, or implicit diffusion
- [Session state 2026-09-10 05:30](session-state-2026-09-10-0300.md) — superseded by the 09:10 note; keeps the FL/V4/V5/V6/B-series detail
- [Session state 2026-09-10 09:10](session-state-2026-09-10-0910.md) — superseded by the 22:00 note: WB-walk overflow + RT neighbour-Planck bugs FIXED and gated (build_guard); lidded run needs no clamp (B13); V9f/V10 open-top tests running to ~12:40; commit plan; user decides prod11 revival / open-top config
- [Session state 2026-09-10 22:00](session-state-2026-09-10-2200.md) — superseded by the 09-11 03:40 note: implicit radial diffusion + angular cap gated; I2 alive past both deaths; I3 (clamp off) and I4 (400 K active medium, from scratch) launched with watcher agents; tree uncommitted (~29 files)
- [Session state 2026-09-10 early](session-state-2026-09-10-early.md) — superseded by the 03:00 note; keeps the FL/V4/catch job list and the agent transcript paths
- [Session state 2026-09-11 03:00](session-state-2026-09-11-0340.md) — superseded by the 12:50 note: implicit radial diffusion + angular cap gated, V8f dropped (Parker wind), I4 active medium
- [SESSION STATE 2026-09-11 12:50](session-state-2026-09-11-0800.md) — START HERE: a10e367d pushed; FOFC 0+1 committed c0163568 (not pushed); I8 done 9e5; CHECK C1 first; uncommitted RT guard (keep) + vacuous dfloor_keep_temperature (drop); decision: top-cell negative density -> FOFC step 3
- [SESSION STATE 2026-09-12 01:00](session-state-2026-09-12-0100.md) — START HERE: FOFC committed locally + viper merged 7cee2eb7 (NOT pushed); gates pass except dhj ck 1-ulp from the RT guard e1db81d8; RG_fofc running, no events at 1.5e5

## Superseded state paragraphs

---

**PREVIOUS STATE (superseded by the 2026-09-13 handover, docs/handover/HANDOVER-2026-09-13.md):** HEAD = fork = fb69c554 on `polar-average-perf`
(bc972cc8 ME fix + MHD implicit conduction; 1bd31298/27ca5b13 cs seam packs bitwise; 77de5618 radimpx1 split, 1 ulp on HIP;
fb69c554 HIP conventions note for orion). cs prod restart 62.0 -> 40.5 ms/cycle; cs grid gap vs sp CLOSED (37.2 vs 38.5).
RUNNING on apu: cs_mhd_prod3 (11611223 + links 11611229/30/31, rot 6 clean, binary swapped to 27ca5b13 for links 2+) and
sp_mhd_prod3 (11569325, rot ~120). s01 AND s05 switch arms all reached 8 rot incl. controls: switches harmless, rescue
UNPROVEN (controls no longer die). NEXT: read cs_mhd_prod3 (drained columns, rt_eiclamp, 1-ME vs sp) at rot ~20-42;
rt_chain_ck is 34% of cs GPU time = next speed lever; the handover list (dt_min 1e-3, density floor vs WB). Details in
[[inflight-2026-09-09-viper]] bottom. Rules: save tokens, delegate to Opus, apudev one-time only, never write in run/,
sbatch from the main session is permitted.

**ORION STATE 2026-09-12 01:00:** the orion snapshot's own bottom paragraph is still dated 09-07 and reads:
"**CURRENT STATE, 2026-09-07 (~09:45), HANDOVER TO ORION.** HEAD 7c652768 on `polar-average-perf`, PUSHED
to the fork (jing-ze-ma/athenak). Viper is in MAINTENANCE 09-07 12:00 -> 09-12 12:00; the next session runs
on ORION, which cannot see viper's bench/ or scratch. Everything needed is in git:
docs/handover/HANDOVER-2026-09-07.md (READ FIRST), docs/handover/scripts/ (the analysis scripts),
docs/handover/claude-memory-2026-09-07/ (a copy of this memory directory as of the handover)."
The live orion entry point is [[session-state-2026-09-12-0100]]: FOFC committed locally + viper merged
7cee2eb7 (NOT pushed); gates pass except dhj ck 1-ulp from the RT guard e1db81d8; RG_fofc running, no
events at 1.5e5. Standing orion rules: [[delegate-heavy-work-to-opus]], [[save-tokens-everywhere]],
[[viper-hip-code-conventions]], [[run-directory-untouchable]], [[test-output-location]].
