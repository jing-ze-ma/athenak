# NOTE 2026-11-24w orion -> Caltech: RSG: check L and R continuity across the hand-over (USER)

USER: check the continuity of the RSG run's luminosity and radius, and fix it if wrong.

The RSG restart from rsg.01707 switches the radiation scheme (old implicit -> clean staggered VET), so L and R can jump. When you start RSG:
1. Write `user.hst` (and `hydro.hst`) EVERY cycle for the first 500 cycles. Add the output dt override on the command line if the deck's hst cadence is coarser.
2. Post the hst + run.log / log.0 after 100 and after 500 cycles to caltech_results/RSG_RC3_1124/.
3. If cheap, a vet + hydro dump at cycle 0 and ~100, so we can measure R_ph (tau = 1 / T_eff radius) directly.
4. Watch the first 100 cycles: if L_top/L_in or the photospheric radius jumps by more than a few % vs rsg.01707's values, or CW / FATAL, pause RSG and post. Our RSGCONT agent compares with the old run's history (on viper) and fixes the cause.

-- orion (coordinator)
