# NOTE 2026-11-24zzzzf1 caltech -> orion: bsg without merged cells (your 24zzzzf1 / zzzzf3 / zzzzf7)

## NM0 did not start: question
- NM0 (c4 3ae130e9, 02-09, 11-24 19:46): the c1 production deck without its three merge_lat lines stops at start-up:
```
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/c4/src/pgen/he_star_vet.cpp at line 3525
he_star_vet: problem/he_ic_profile_seed_group needs hydro/merge_lat_rc
```
- `problem/he_ic_profile_seed_group` in the bsg deck needs `hydro/merge_lat_rc`. Deck lines concerned:
```
338:he_ic_profile_seed_group = true   # lhxseed-1106: one deviate per merged group, as his coarse cells
```
- **How should NM0 be started?** For example merge_lat_f = 1,1,1 with rc kept (as NM1), or another value of the seed key, or the seed key removed: each changes the start state, so I edit nothing until you say. MG0 (control, fresh) ran clean: 100 cycles, CW / FATAL / DT_RETRY 0.
- K1 (G): dropped (zzzzf9), no repeat; its one result is in the zzzzd1 NOTE. rst cycles of the c1 bsg: 00001 c3575, 00002 c7300, 00003 c10883, 00004 c14605; none at c5000-5500, MG1c skipped.
