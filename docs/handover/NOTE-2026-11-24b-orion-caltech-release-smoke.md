# NOTE 2026-11-24b orion -> Caltech: release-branch GPU smoke (B = IA12.2 vs R = cleaned release), after WD

The release branch is the cleaned code that others will use (m1 -> vet rename, staggered only). Before we publish it, it must be bitwise to IA12.2 on H200. This is a quick B vs R check, run AFTER the WD rerun (NOTE 24a).

## Binary R
- Branch `caltech-release-38f3d3c95` (src = cc7a47d4a), CUDA 12.9 + g++ 13 (never 12.2).
- Build with `PROBLEM=he_star_m1`: this tests the cmake alias. B is your existing athena_ia122_465aff485.

## Arms (B vs R, identical command lines with the OLD key names, hst output only, no rst kept)
1. bsg: `ia122_decks/tile/IA122_bsg.athinput`, fresh, `time/nlim=30`, 4 GPUs.
2. m35mms: `ia122_decks/tile/IA122_m35mms.athinput` from your latest IA12.2 m35mms rst, +20 cycles, 4 GPUs.
3. t9l62: `IA122_t9l62_restart.args` from rst 00012, +10 cycles, 4 GPUs. This also checks that R reads an rst that stores `<rad_m1>`.
4. (only if you have old-prod HG / RSG rsts; otherwise skip.)

## Check per arm
- rc 0, 0 FATAL, the same CLEAN WARNING count as B;
- the non-# rows of every .hst identical (cmp / md5) B vs R;
- R's log has no "Block name 'rad_m1' ... not found" and no "never read" key that B reads. R may list cs_cross_mono / implicit_pos_floor_solve / vet_gd_reach as "no effect"; that's expected.
- Optional: the docs' STARS.md sp example on 4 GPUs, `time/nlim=400`: rc 0, 0 CW, `clean_init steady pass 1` = 9.025e-09.

Post one table (arm, B md5, R md5, equal?).

-- orion (coordinator)
