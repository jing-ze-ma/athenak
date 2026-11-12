# NOTE 2026-11-13 Orion -> viper: consfresh. Part (b) d6fast3 fresh-start gate: **PASS** (part (a) audit follows)

## Builds (git archive + kokkos from the Orion src_vfix copy; patches applied with git apply, both clean)
| tree | commit | patch | athena md5 |
|---|---|---|---|
| f | fbb5c663825e3861140c34430dc5f804b6d45200 | - | dbbaf51588ab44268ab0aa8e2987cd1b |
| p | d3b3f2e5531a915c69b8ede3347fd1aca8951860 | - | 72ccef10641b7ad34e1b355d246f7fb3 |
| fa | fbb5c663 | consaudit_on_fbb5c663.patch | 8ddfeb83f3f4bfe59fd1e6c309d64019 |
| pa | d3b3f2e5 | consaudit_on_d3b3f2e55.patch | cb82bd1f82190da10a8896e1ed81d781 |

Build: mpicxx gcc/13 + openmpi/4.1.8, MPI ON, PROBLEM=he_star_m1, Release, -ffp-contract=off. `strings fa | grep -c dbg_cons_audit` = 1.
Tables: hgpos-1109/files, md5 ic f9678c2a, planck 23ee543f, rosseland e2868b75 (= viper). 24 ranks per run, `mpirun --bind-to none`, binding off, OMP 1.

**Deviation (needed):** the $AM and $DG keys cannot go on the command line, because none of them is in the decks. AthenaK FATALs with
"Parameter 'am_cons_diag' in block 'hydro' on command line not found" (and the same for dbg_energy_audit). I inserted them into each run's copy of
the deck directly after the <hydro> / <rad_m1> header (addkeys.py). time/nlim stays on the command line (nlim is in the deck). Smoke (3 cycles,
4 runs) then gave rc 0 and 0 FATAL, with AMCS_CONS / CAUDITH lines printing.
G and N ran on orion060 (4 x 24 ranks), A on orion101 (2 x 24), p.exclusive.

## (b) Gate, 40 cycles from the IC (plain binaries + AM print keys)
| run | rc | FATAL | NOT ADMISSIBLE | NOT POSITIVE | stage / BE / fallbacks | Picard solves, mean / max / NC | wall s |
|---|---|---|---|---|---|---|---|
| G_D6 | 0 | 0 | 0 | 0 | 39 / 1 / 0 | 7.900000e+01, 3.278481e+00 / 5.000000e+00 / 0.000000e+00 | 288 |
| G_PROD | 0 | 0 | 0 | 0 | 39 / 1 / 0 | 8.300000e+01, 3.265060e+00 / 5.000000e+00 / 0.000000e+00 | 272 |

- **D6 gate PASS** (rc 0, FATAL 0, NOT ADMISSIBLE 0, NOT POSITIVE 0; 39 stage steps, 1 BE, 0 stage fallbacks).
- **face_fct (hesdirk2 summary, D6):** calls 78, **limited cells 50881**, max rounds 8, withheld E x vol 1.351e40, left (Q<=0) 0.
  face_mmat: checked passes 259, limited cells 0. PROD has no face_fct.

### amlz_analyze.py G_D6 G_PROD
```
G_D6 AMH stages 80 sum_dT_z 2.371e+38 max|resid_z| 2.535e+29 Lz 1.077552e+29 -> -4.619820e+38
G_D6 AMCS solves 79 max|resid_z| 1.958e+26 max|resid_xy| 2.446e+39 mean scale_z 8.771e+40
G_D6 AMCS_DIAG sums |dLr_z| |dLg_z| |dtTpl_z| |dtTrad_z| |del_z| |del_z|/(rhoVR) a_f |Lg_z|: 6.297e+40 4.859e+40 8.929e+40 6.486e+40 4.770e+39 2.791e+01 9.948e+38 5.092e+41
G_PROD AMH stages 80 sum_dT_z 1.658e+38 max|resid_z| 2.080e+29 Lz 8.254908e+28 -> -2.531348e+38
G_PROD AMCS solves 83 max|resid_z| 5.077e+25 max|resid_xy| 2.446e+39 mean scale_z 1.126e+41
G_PROD AMCS_DIAG sums |dLr_z| |dLg_z| |dtTpl_z| |dtTrad_z| |del_z| |del_z|/(rhoVR) a_f |Lg_z|: 6.091e+40 4.872e+40 9.081e+40 6.723e+40 4.789e+39 2.796e+01 -3.050e+38 5.098e+41
```

- **gas L_z after 40 cycles: D6 -4.62e38, PROD -2.53e38** (viper GPU PROD at 30 cycles: -4.5e38). **AMCS_DIAG |del_z|: D6 4.77e39, PROD 4.79e39**
  (viper GPU PROD 30 cycles: 2.2e39; 40 cycles here). D6 has the production bookkeeping: |del_z| equal to PROD within 0.4 %, |L_z| within 2x of PROD
  and of the same order as viper's production reference. AMH max |resid_z| 2.5e29 (D6) vs 2.1e29 (PROD); AMCS max |resid_z| 2.0e26 vs 5.1e25.
- **Neutrality, part 1 (G vs N, 40 rows):** the plain build with the AM print keys (G) = the diagnostic build with both keys off (N), md5-identical:
  D6 hydro 0b2bef81 = 0b2bef81, user e2be5729 = e2be5729; PROD hydro 48502d66 = 48502d66, user e458e658 = e458e658.
  The A (keys on) comparison follows with part (a).
