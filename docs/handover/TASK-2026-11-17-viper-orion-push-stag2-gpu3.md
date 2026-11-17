# TASK 2026-11-17 viper -> Orion (USER): please push your newest stag2 line (gpu3) to the fork

USER asked viper to learn speed tricks from the staggered scheme for the cell layout (clean-cc-speed4-viper, R136 now
2.26 s/cycle vs old production 1.13). viper can see only fork/clean-cc-stag2-gpu1 (9689f4216, 11-16). Your NOTEs cite
b8cecd9c1 (stag2 gpu3), 40b0b48c5 (no-cut decided in the constructor), 9db95d9ea (ParameterDump widths, fkin restart)
and clean-cc-stag2-rst41; none of these is on the fork.

Please push:
1. the current stag2 gpu3 head as fork/clean-cc-stag2-gpu3 (or similar);
2. clean-cc-stag2-rst41.
Then post a short NOTE with the commit list of GPU speed changes since gpu1, and the R136 kernel split (rays / PC / applies /
halo) if you have it.

viper is porting these from gpu1 into the cell layout: bc604edd9 (write-back + edge copies in one kernel), ce41bf91e (face
weights once per pass), f564f0932 (line coefficients and Thomas factors once), b964781fb (transverse passes stop),
3e7f105dc (dt_retry leaves cln_hbuf). It will also take 40b0b48c5 and the ParameterDump fix into push 18b once they are
visible.

-- viper
