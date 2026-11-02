# NOTE 2026-11-02 viper -> Caltech: roche_wb_tan default on for the cs accretor merged; fin51-r-fix = c38bd1981

This note supersedes the accwbtdef-1102 part of ee4d29826. The fdb167169 merge on viper was never pushed.

## New fin51-r-fix = c38bd1981
- c38bd1981 = merge --no-ff of accwbtdef2-1102 (4c5f540f4) into 24d044fc1. It contains:
  - ab4370a8d: problem/roche_wb_tan default logic.
  - f705954be: your 2a71977a3, canon decks set the key explicitly.
  - 4c5f540f4: your 54e74ee0b, the finalize fix (RyPerFinal frees wbt_ga2_/wbt_ga3_).
- **It merges cleanly into prod-1102:** `git merge-tree --write-tree fork/prod-1102 c38bd1981` reports no conflict.
- **Default rule when the key is absent:**
  - true when all of these hold: cs, fresh start, phys_3d, no spin ramp, etotgrav, (no env_ic_vz or thermo general), single level;
  - otherwise false. The sp accretor and restarts whose input lacks the key stay false.
- A named key is used as given.
- One root log line gives the choice, and the restart input records it.
- **roche_wb_tan is default on for the cs accretor (fresh start). Please take it into prod-1102.**

## MI300A gates
ry_per_accretor HIP builds, apudev, 2 ranks, 10 cycles. Deck: cs16 (canon_cs96_static_s1 at 16x16 per panel) with column_mesakap. Base is 24d044fc1. Smoke job 12608230, gates job 12608231.

| gate | result |
| --- | --- |
| (1) new: cs16 key absent vs key = true | rst data IDENTICAL; hydro/trc/user/wbt hst IDENTICAL |
| (2) cs16 key absent: base vs new | DIFFER as required; log says "true (cs default, fresh start)" |
| (3) sp run12 deck (nx3 128), key absent: base vs new | rst IDENTICAL, hst IDENTICAL |
| (4) cs16 restart, key line stripped, +10 cycles: base vs new | rst IDENTICAL, hst IDENTICAL; log says "false (restart input without the key)" |
| (5) dhj d_sts 20 cycles | rst and mhd.hst IDENTICAL; the dhj binary md5 equals 34c7b76e's |
| key-on runs | rc 0, no "deallocated after Kokkos::finalize" |
