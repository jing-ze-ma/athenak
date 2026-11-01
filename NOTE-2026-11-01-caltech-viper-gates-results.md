# NOTE 2026-11-01 Caltech -> viper (USER): bitwise gate results (G3 TRC PASS; G4 HALO bitwise PASS but no saving: fallback on 99.9 % of super-stages)

Answers TASK 4bdc9e76 (gates moved). G1 POL and G2 CSV were dropped here because viper passed them (c3d73bb2). G4 HALO (TASK 6e5e226a part B) is at the end.

## Binaries (H200, fmad off, builds_rstfixx/bin)
| tag | commit | binary | md5 |
|---|---|---|---|
| BASE | fin51-r-fix 65658289 | athena_gpu_dhj_65658289_fmadoff | e33f2e6abaf57573514ab5c7bde5c1b5 |
| TRC | tracers-defaults-1101 8ccae663 | athena_gpu_dhj_trcdef_8ccae663_fmadoff | de6bf3f204de89d84ef0ebce296d707b |

All runs: 1 rank, 1 H200 (hpc-sm-02-03), BVALS_CC_SPLIT=1, output1/dt=1 (hst every cycle), final rst. Comparison: rst bytes after `<par_end>`; dhj.mhd.hst as text; dhj.user.hst per-column max relative difference.

## G3 TRC: PASS
| sub-gate | start | cycles | rst payload | dhj.mhd.hst | dhj.user.hst |
|---|---|---|---|---|---|
| (a1) old rst WITHOUT tracer keys, BASE vs TRC | t0a_ref dhj.00002.rst (r1m_prod_nosplit deck, written by b935e145; no xe/floor_tracer keys; t = 2000 s, cycle 186) | 200 | IDENTICAL | identical (181 rows) | only Efloor/Mfloor/Efloor_rt (cols 12-14) differ, <= 2.3e-14 |
| (a2) developed 10 G rst (keys explicitly false), BASE vs TRC | keep/off_0p50.rst (T1 off arm, t = 0.5 rot, cycle 8774) | 200 | IDENTICAL | identical (200 rows) | cols 12-14 only, <= 1.2e-13 |
| (b) fresh start, TRC: no keys vs explicit xe_tracer = true, floor_tracer = true, valfven_noinfall_fac = 4 | r1m_trc deck, t = 0 | 150 | IDENTICAL | identical (151 rows) | cols 12-14 only, <= 1.9e-14 (25 columns in both: the tracer columns are present) |

- Log of (a1) with TRC: `<mhd>/floor_tracer = off (restart input has no key: the old default)` and the same for xe_tracer.
- Log of (b) without keys: `valfven_noinfall_fac = 4`, `floor_tracer = true`, `xe_tracer = true` reported as defaults.
- The user.hst differences are the known GPU-atomic floor sums (columns 12-14 only), far below 1e-12.
- Note: the (a1) rst has an explicit valfven_noinfall_fac = 0.0 in its deck, so the fac default is exercised by (b) only.

Run dirs: /resnick/groups/carnegie_poc/jingze/trc_1101/runs/g3a1_*, g3a2_*, g3b_*. Decks for (b): /resnick/groups/carnegie_poc/jingze/gates_1101/g3b_{nokeys,explicit}.athinput.

## G4 HALO (rkghalo-test-1101 @ 0594b145): bitwise PASS, but no saving on r1m: the fallback fires on ~every super-stage
Binary athena_gpu_dhj_rkghalo_0594b145_fmadoff, md5 0ac9ee80957f38ea1a69b1d3189092f2. Start keep/off_0p50.rst (developed 10 G r1m, cycle 8774), production keys except `mhd/rkg_lean=false`, 1 H200 (02-03 GPU 3, alone for all G4 runs). RKG 16 stages/cycle in this state.

| test | result |
|---|---|
| 1. smoke, 10 cycles, key on | runs clean, no FATAL; `rkg_halo_ien_only: 144 super-stages >= 2 fell back` (of 150) |
| 2. bitwise, 200 cycles, off vs on, hst every cycle | rst payload IDENTICAL; dhj.mhd.hst identical (200 rows); dhj.user.hst: only Efloor/Mfloor/Efloor_rt differ, <= 2.1e-13 -> PASS |
| 3. fallback count | key on: 2994 fallbacks in 200 cycles x 15 super-stages >= 2 = 3000 (99.8 %); in the 300-cycle cost runs 4494 / 4500 (99.9 %) |
| 4. cost, 3 x 300 cycles interleaved off/on (s/cycle over cycles 8800-9000 from the run.log stamps) | off 0.2643 / 0.2631 / 0.2625 (mean 0.2633); on 0.2682 / 0.2680 / 0.2676 (mean 0.2679): **on is +1.7 % slower** |

- Verdict: correct (bitwise), but useless on the r1m state, because some ConToPrim (the floors) writes a non-energy slot on almost every RKG stage. The Allreduce check plus the full exchange costs a little more than the plain path.
- To get the 43 ms/cycle, the fallback test would have to ignore the floor-only writes, or the floor ConToPrim would have to be excluded from the RKG sub-stages. Either changes the answer, so it is no longer exact.
- Run dirs: /resnick/groups/carnegie_poc/jingze/trc_1101/runs/g4_* (script gates_1101/g4.sh).
