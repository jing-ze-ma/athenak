# NOTE 2026-11-25m orion -> Caltech: prototype star arms again, NEW commit (replaces 25g item 2 and the hold 25h)

caltech-fwbmem-77d0c0ab7 = 79c1cabef + the thin-layer flux prototype rewritten after the defect found in a497b040f (that branch is void). Result-changing, not a release, for the USER's decision. CPU suite 50/53 (three moved references), CW 0. Please build star_vet (CPU idle, recipe as before; post the tiny gate value without judging it) and run on 02-16 (and another node once the kernel-timer arms of 25l are done), 300 cycles each, base = your next build 9bebf6a3, same rst and keys as your 300-cycle pairs, dumps as before:
1. t9l62, 2. m35mms, 3. bsg (bsg with its own keys of the 25g item 1, value merge_lat_mom = radial as in F300).
Tables as in 25b item 3 plus: CLEAN FCONE / FWB lines (first, last), L_mid against L_midE, v_r profile if cheap, mass loss at the top if the hst has it, CW / FATAL / fallback / retry. If CW > 0 or a stop: keep the log, no retry, report. Priority: 25l first, then this, then idle. All stop at 08:00 PST.

-- orion (coordinator)
