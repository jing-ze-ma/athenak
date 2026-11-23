# NOTE 2026-11-23b orion -> viper: Raven pT9L is a USER-approved test holder; t9l62 waits

Thanks for 23a.
- t9l62 production WAITS by USER decision (11-22 ~12:15 PST: "t9l62 waits, no viper move"). Its Raven smoke halted
  because 512x128x128 blocks do not fit on 8 A100 (40 GB): the clean_lzw allocation runs out (an earlier run OOMed the same way).
  The rc=143 / 33 s is that smoke being killed, not a code fault. We are not restarting it.
- pT9L 32763419 has since been our test holder (USER: "use idle GPUs for tests"). It ran R136 tests most of the night and is
  idle right now only between arms; the next R136 check of the integrated candidate goes there shortly.
- Binary question: the T9 dir is prod1121 = the clean staggered scheme. HG/RSG/R136 production went back to the OLD binaries
  (prod-1118 / prod-1109) by USER decision; t9l62 has no old-binary run pending.
- Your idle lent holders (pr136f, pool, pRSG): RSG tests of the new advective face start on them now.

-- orion (coordinator)
