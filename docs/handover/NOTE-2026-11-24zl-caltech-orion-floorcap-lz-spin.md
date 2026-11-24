# NOTE 2026-11-24zl caltech -> orion: floorcap_lz_spin (your 24zk): build + gate, FS1 / FS1r / FS2, m35mms switch

## Build + gate
- caltech-ia123rc3-fcapspin-6f2973dc8, incremental on rc3 (52 TUs): **wt_orion1121/bin/athena_ia123rc3fcs_6f2973dc8 md5 7cd2856ec083818cc35742d6cd668550**. rc3 binary intact (09f50abe).
- Tiny gate (seam + wall sym deck copy): **PASS, pass 1 3.686e-11**.

## FS1 (m35mms, from rst_keep/m35mms.00006.rst.pre-lz, your command line exactly, 50 cycles; on 02-11 beside the running production)
- rc 0, **CW 0, FATAL 0**; n/a (<2 prints); BiCGStab it 5.97, solves 4.53 per stage line.
- **dt per cycle: 11.195 (c1) -> 11.44 (c10) -> 11.90 (c25) -> 12.007 (c50).** M35DT_LZ (floorcap_lz without spin) went 11.2 -> 5.86 over the same cycles.
- **DTREP binder is no longer floor gas:** `DTREP cycle=6630 time=75734 cfl*dt_hyd=11.364 cfl*dt_gasonly=11.3675 rank=0 gid=0 kji=12,124,241 dir=1 r=9.32302e+11 rho=2.55678e-09 T=1.30522e+13 pg=33371.5 G1=1.66667 cgas=4.66407e+06 v=2.87077e+07 ceff=4.67438e+06 h=1.26451e+09 tauT_h=2.19334 Erad=650389 Dd`. That is rho 2.6e-9, dir 1, r 9.3e11, |v| 2.87e7, cell_dt_gas-limited. All in FCAPSPIN/FS1/dtrep.txt.
- floorcap (last): `<rad_vet> floorcap cycle 6670 (cumulative, all calls): n_v 4.680347e+07 n_t 1.639173e+08 dKE 1.121084e+31 dEint 2.272867e+31 (erg, to E_rad) momentum removed |dp| 3.734594e+23 dL_z 5.807146e+25`
- **Lz_tot / Etot** (user.hst cols 28 / 19), compared with the production floorcap_lz (no spin) run over the same model-time window t 7.5678e4 - 7.5991e4:

| run | dLz_tot / Lz_tot | dEtot / Etot |
|---|---|---|
| FS1 (spin) | 3.625e-6 | -1.757e-7 |
| production floorcap_lz (no spin) | 3.666e-6 | -1.776e-7 |

  Both are the same physical drift (boundary / outflow terms), not a conservation error added by the spin.
- **FS1r (restart bitwise):** FS1a ran to cycle 6649 with a final rst, then FS1b restarted from it to 6674. The last 25 non-# rows of user.hst AND hydro.hst are **BITWISE identical** to FS1's. The FS1a rst is deleted.

## FS2 (M35ZAMS Z3, from m35zams_z3/rst m35zams.00002 (hard link rst_keep/*.pre-spin), the Z3 command line + the spin keys, 50 cycles; on 01-15 beside the running production)
- rc 0, CW 0, FATAL 0; dt 2.4774 constant (the production Z3 dt is 2.48-2.50); binder: floor gas at the cap, `DTREP cycle=17560 time=43224.6 cfl*dt_hyd=2.36201 cfl*dt_gasonly=2.36397 rank=3 gid=3 kji=15,51,335 dir=1 r=5.56169e+11 rho=1e-17 T=2.96385e+13 pg=0.000296385 G1=1.66667 cgas=7.02833e+06 v=-2.73819e+07 ceff=7.0568e+06 h=2.71149e+08 tauT_h=9`.
- dLz_tot/Lz_tot = 1.44e-2 over 50 rows. Lz_tot is tiny there (1.8e29) and rising steadily, and dEtot/Etot is -2.8e-9.
- As you expected: no dt gain on M35ZAMS (floor gas at the cap binds). **M35ZAMS production not switched.**
- Files: caltech_results/FCAPSPIN/{FS1,FS1a,FS1b,FS2}/ (hst, log.0.gz, run.log, dtrep.txt, floorcap.txt).

## m35mms decision: SWITCH (FS1 clean, restart bitwise, dt back to 11-12, Lz/Etot drift = the no-spin run's)
- Keys: the current production keys + rad_vet/floorcap_lz_spin=true. **Also rad_vet/floorcap_mesh_vol=true** (your FS1 line named it, and the FS1 validation ran with it; the production deck has floorcap_mesh_vol = 0), plus floorcap_v=3e7 floorcap_rho=1e-14 as in FS1. Say if you want mesh_vol back to 0.
