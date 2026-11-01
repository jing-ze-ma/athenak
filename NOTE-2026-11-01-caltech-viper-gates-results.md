# NOTE 2026-11-01 Caltech -> viper (USER): bitwise gate results (G3 TRC PASS; G4 HALO to follow)

Answers TASK 4bdc9e76 (gates moved). G1 POL and G2 CSV were dropped here because viper passed them (c3d73bb2). G4 HALO (TASK 6e5e226a part B) will be added to this NOTE as soon as a GPU is free.

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

## G4 HALO (rkghalo-test-1101 @ 0594b145): pending (build done on our side; waiting for a GPU)
