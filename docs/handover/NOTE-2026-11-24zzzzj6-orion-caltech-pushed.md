# NOTE 2026-11-24zzzzj6 orion -> Caltech: release pushed (ccab558ee); L and R of the restart tests; branch clean-up

1. The release is on GitHub: `vet-clean-release` = ccab558ee (its src is b09a7b3c1 plus comment lines). Thank you for the suites. Continue the t9l62 / t9mlt chains.
2. USER: the main thing to check for a restart is whether **luminosity and radius are continuous and smooth**. For each star of zzzzj4 (R50 vs C50 from the same rst, and the production rows before / after): L_top, **L_bot**, L_mid, L_int at the first step, steps 1-10 and 50, and the same for the run before the restart file; jump at the first step and slope. Radius: from dumps at cycle 1, 10, 50 (lateral-mean tau = 2/3 radius by the run's opacity, T_eff); if there are no dumps say so. Flag any jump above the run's own noise. (On RSG L_top is continuous to 1e-3 but L_bot jumps +4.6 % at once and +11 % by cycle 5: old-scheme start; look for the same in yours.)
3. When your tests no longer need the `caltech-*` branches on GitHub, say so: they will be deleted (USER).

-- orion (coordinator)
