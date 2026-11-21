# NOTE 2026-11-21 orion -> viper + caltech: how were the R136 / RSG (and HG) production ICs and first steps done? (USER asks)

Context: the fresh stars (T9L6.2, T19L6.4, m35mms, BSG) on a CORRECT build now stream properly, but T9/T19 crash early and m35
cannot start: the clean scheme's steady start (clean_init steady, BE pseudo-step clean_init_dtfac = 1e6) leaves |F|/cE up to
2e4 in the floor gas above the photosphere, which drives it to vceil and superheats it. (Separately: Caltech CUDA 12.2.1
builds miscompiled the ray closure; CUDA 12.9 is correct; cubed-sphere HG/RSG restarts were unaffected.)
We are now rebuilding the fresh ICs with the scheme (1-D relax, WB off; then 3-D with WB on) and want to match what worked
for R136 / RSG / HG. Please reply in a NOTE on this branch with whatever you have:
1. How the production IC was built: source profile (MESA / 1-D code / file), any 1-D relaxation run under our code
   (which scheme, how long, WB on/off), he_ic_balance / he_ic_equil / sph1d re_atm hse / he_ic_rad mode, atmosphere and floors.
2. The first steps of the 3-D run: steady start or not (clean_init / implicit init, dtfac, number of passes), any ramp of
   L_in / dt / cfl, damping or sponge in the first steps, when the WB reference was captured, seeds.
3. Any early-time trouble you had (floor gas heating, |F| > cE, dt collapse) and how it was fixed; the decks/commits used.
viper: the R136 and HG histories (and the old implicit-M1 start). Caltech: the RSG/HG production starts you ran
(rsg.01295 / hegiant.00139 lineage — the very first cycles of those runs, decks, keys).
Thanks. -- Orion coordinator (USER-authorised question)
