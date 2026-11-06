# NOTE 2026-11-06 viper -> Caltech, orion (USER): reference setup is what Linhao's code and deck say

**USER decision: use Linhao's own code and deck as the reference, including GM = 260.795 code (19.78 Msun, read from `<problem>/GM` in his deck), not 20 Msun.**

## The reference values
From `athenapp-bsg/linhao_input/athinput.bsg_star.rusty` and his `bsg_star.cpp`:
- **Gravity:** GM = 260.795 code. No self-gravity.
- **Time step:** cfl_number 0.3, **cfl_rad 0.025**, so dt ~5 s (`new_blockdt.cpp:200`).
- **Units:** density_unit 6.15245e-5, T_unit 1.07887e6, length_unit 6.957e10, molecular_weight 0.62.
- **Radiation:** nmu 4, angle_flag 1, nzeta 5, npsi 6, taucell 4, error_limit 1e-6, nlimit default 100.
- **Floors:** dfloor 1e-9, pfloor 1e-12, tfloor 1e-4 (code).
- **Seed as coded:** amp 0.05 where rho > 0.01 code (6.15e-7 g/cc).
- **Grid:** 192 x 128 x 128, 16^3 blocks, SMR level 1 at 38-55 (2560 blocks), x2 and x3 periodic.
- **Tables:** his original `BSG/log_*.txt`.

The only non-physics differences in our runs: CRs stripped from the deck, and the 4-line `rad_source` user-output edit in the pgen.

## Orion
No change is needed. Your `prod_seeded` copy (job 227350) already uses this deck, so keep it as queued. Raven's copy (31975424) is running with the same setup.

## Caltech
Please use these values when comparing with Athena++, in particular GM and cfl_rad in your 1-D and any AthenaK reference arm. Our earlier arms at GM = 20 Msun and cfl_rad = 1 (A_af1, A_tc4, `af1_cflrad`, audit, ledger) are not his setup. Use `his1d` as the 1-D reference (`athenapp-bsg/cflrad_1106/`).

## Status of the caution in NOTE-2026-11-06h
Whether the published run used this same GM remains unconfirmed, because his Fig. 2 KE near 1 d (~7e40) sits nearer the 20 Msun arm (5.2e40) than his1d (5.2e39). We proceed with his code and deck regardless; the 3-D run will show whether they reproduce the figure.
