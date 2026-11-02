# NOTE 2026-11-02 viper -> Caltech: csaccmerge-test-1102 merged into fin51-r-fix

## 1. csacc merge: fin51-r-fix = 24d044fc1 (pushed)
- 24d044fc1 = merge --no-ff of 34c7b76e8 into 7f0e3dd90. The tree is identical to 34c7b76e8.
- MI300A gates: dhj_gpu72 build, apudev, keys absent, base 7f0e3dd9 vs test 34c7b76e8. Job 12606273 (smoke 12606161).

| arm | deck | cycles | ranks | rst data after par_end | mhd/hydro hst |
| --- | --- | --- | --- | --- | --- |
| dhj MHD | d_sts (r1m-style w121 MHD, STS) | 40 | 2 | IDENTICAL | IDENTICAL |
| dhj MHD | d_nosts | 40 | 2 | IDENTICAL | IDENTICAL |
| dhj MHD | d_sts | 40 | 1 | IDENTICAL | IDENTICAL |
| W121 dhj.00641.rst | +20 | 20 | 2 | IDENTICAL | IDENTICAL |
| dhj hydro cs | r1m_viper_hydsc, am_hst_split=false | 20 | 2 | IDENTICAL | IDENTICAL (user.hst too) |
| base vs base repeat | d_sts and hydro cs | 40 / 20 | 2 | IDENTICAL | IDENTICAL |

- The MHD user.hst Efloor/Mfloor/Efloor_rt columns differ at <= 8e-14 rel. They differ just as much base vs base, so the cause is a non-deterministic reduction, not the merge.

## 2. accwbtdef-1102 (9fb488887): roche_wb_tan default on for the cs accretor (fresh start)
- **Status:** gated, and the branch is pushed to fork. The merge into fin51-r-fix is prepared on viper (fdb167169) but **not pushed yet**: it waits for the user's OK. Until then fin51-r-fix stays 24d044fc1.
- **What it changes:**
  - If problem/roche_wb_tan is absent, it defaults to true when all of these hold: cs, fresh start, phys_3d, no spin ramp, etotgrav, (no env_ic_vz or thermo general), single level. Otherwise it defaults to false.
  - The sp accretor and restarts whose input lacks the key keep false.
  - One root line in the log names the choice.
  - Caltech 2a71977a3 is cherry-picked (canon decks set the key explicitly).
- **Bug fix:** wbt_ga2_/wbt_ga3_ are now freed in RyPerFinal. Before, every key-on cs run aborted after Kokkos::finalize (rc 134, after the last output). This affects runs on csacc-1031 / 24d044fc1 with the key on.
- **Gates:** ry_per_accretor HIP build, apudev, 2 ranks, 10 cycles. Deck: cs16 (canon_cs96_static_s1 at 16x16 per panel) with column_mesakap. Jobs 12607289 and 12607656 (smoke 12607149).

| gate | result |
| --- | --- |
| (1) new, cs16 key absent vs key = true explicit | rst data IDENTICAL, all hst IDENTICAL |
| (2) cs16 key absent, base 24d044fc1 vs new | DIFFER as required; log says "true (cs default, fresh start)" |
| (3) sp run12 deck (nx3 128), key absent, base vs new | rst IDENTICAL, hst IDENTICAL |
| (4) cs16 restart, key line stripped from the rst header, +10 cycles, base vs new | rst IDENTICAL, hst IDENTICAL; log says "false (restart input without the key)" |
| (5) dhj d_sts 20 cycles | rst and mhd.hst IDENTICAL; the dhj binary md5 equals 34c7b76e's |

- **Request:** roche_wb_tan is default on for the cs accretor (fresh start). Please take it into prod-1102 once fin51-r-fix carries it.
