# NOTE 2026-11-24zzzzc5 orion -> Caltech: t9l62, is the outburst one-time? (read-only; USER question)

USER's hypothesis: the fresh t9l62 outburst may be a start-up transient of the radiative-only IC (the IC is Jiang's init_outburst.txt with F_r = L_IC / (4 pi r^2) everywhere and no convective flux), not a property of the star. Test on EXISTING histories; no run, nothing written in production directories. Low priority, after zzzzc4 and behind all production and arms.
1. List every t9l62 run you hold with a user.hst (fresh prod 1120 / 1121b / rc3 / T2 / rc4v, S5x, PL etc.): name, binary, scheme, start (fresh IC or which rst), time range.
2. Per run (if columns exist): L_top / L_in, L_bot / L_in, KE total, M_tot (mass leaving the top: integrated outflow, or M_tot change), max or mean radial velocity if a column exists, dt, versus time, on a 2e4-s grid (text table), and the long fresh run (to the latest t you have, 1.6e6 or more, across restarts where the scheme changed: mark them).
3. Count episodes: an episode = dt < 0.5 x its quiet value (96.85) for more than 10 consecutive cycles, or KE or L_top excursions you would call outbursts; give start and end times, minimum dt, peak L_top and peak KE of each, per run. State that you used this definition; it is mine and arbitrary.
4. Where the fastest outward mass sits at 2e5, 4e5, 4.7e5, 5.7e5, 8e5, 1.2e6 s (if dumps exist): the radius of the max mass flux 4 pi r^2 rho v_r, and whether Mdot decays after the first episode.
5. Do the runs differ in their first episode's start time / strength (same IC, different keys)? A table only; no interpretation needed.

-- orion (coordinator)
