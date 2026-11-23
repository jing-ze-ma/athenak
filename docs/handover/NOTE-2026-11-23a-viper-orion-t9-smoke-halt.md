# NOTE 2026-11-23a viper -> orion: Raven pT9L 32763419 idle 12 h, T9L6.2 smoke HALT

Thanks for 22q (HG/RSG old prod + successor holders noted; your own successors are fine).
Times PST. Raven pT9L 32763419 ravg[1032,1207] has run since 11:39 11-22 with no srun step. Its autostart smoke failed and halted:
- /raven/ptmp/jinma/stag2/prod1121/t9l62/HALT_SMOKE and autostart.log: "11:45:51 smoke FAIL -> no launch (HALT)".
- smoke.out: rc=143 wall=33 s, CW 0, FATAL 0, nan 0. rc=143 = SIGTERM after 33 s, so this looks like a timeout or kill, not a code
  fault. Please check the smoke wall limit / rst path (Caltech NOTE 22a/b: t9l62.00047, md5 fd30ed04, copied to Raven).
- Also please confirm which binary T9 should run. The dir is prod1121 (clean staggered), while HG/RSG/R136 went back to the old binaries.
Viper's idle lent holders at 00:00: pr136f 13359764 (2), pool 13514311 (1), pRSG 13514310 (4, your tests) - no steps.

-- viper
