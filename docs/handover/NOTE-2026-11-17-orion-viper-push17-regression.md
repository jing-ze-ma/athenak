# NOTE 2026-11-17 Orion -> viper (USER): push 17 e5c54306a CPU regression, PRELIMINARY (NOT GREEN yet)

Answers TASK-2026-11-17-viper-orion-push17-regression.md. **Status: not green. A bisect is running;** this file
will be updated with its result.

- e5c54306a: md5 `cf52f1cec16c3441acc78ab5e00782c4` (git archive build, Orion gcc/13 Release, -ffp-contract=off).
- dd468ab99 (push 16): md5 `824b439cbb0bcdf6d241a781357be5d7`.
- Everything ran in `srun --overlap -c 112 --cpu-bind=none` steps of hold 230761 (orion002).

## (A) defaults vs dd468ab99: bitwise everywhere except the cubed-sphere RT feed deck
- **Bitwise** on the reg harness (30 runs: bin data, tab/hst, CLEAN lines). This includes multi-rank, the x1line RT
  deck and WB fv1/fv/lld/cart (rh.user.hst byte-identical).
- **Restarts bitwise** for hcyl, sp3mb, sp3mbmr, thin1, rot1s and csCnr.
- **Not bitwise: RT feed / feedx1 / feedmr** (rst/feed.athinput: cubed sphere, fv1, `clean_face_model = milne`,
  `force_reference = clean`).
  - The **fresh run changes**: 16-cycle maxrel 2.8e-3 grouped; m1_e 1.65e-4, m1_f1 1.9e-3.
  - **Its restart is not bitwise** either: 16 straight vs 8 + restart + 8 differ by m1 7.9e-4 abs and hydro_u
    1.1e-5 abs.
  - The change is reproducible run to run.
  - With `clean_halo_vars = 32` in the deck the results are the same, so the halo default is not the cause.
  - cs/feed1 (the same cubed-sphere feed without milne and force_reference) is bitwise.
- 0 CLEAN WARNING, 0 DT_RETRY, 0 FATAL.

## (B) opt-in arms (keys written into the decks; 16 cycles; maxrel vs the default straight run)
| deck | clean_fv_lag = extrap | float set (fv_float, pc_factor_float, extrap_share, dt_retry_views = clean) | float set + forced retry (dt_retry_force_cycle 4, fac 1) vs float set |
|---|---|---|---|
| sp3 (sp, fv1, force_reference = clean, 2 blocks) | 4.0e-5 | 6.2e-8 | **6.4e-6, not bitwise** |
| thin1 (cubed sphere, fv1) | 2.7e-12 | 3.9e-9 | **5.7e-16, not bitwise** |
| hcyl (Cartesian, clean_rays = fv1; P17 = P16 bitwise) | 6.4e-16 | 1.3e-6 | **bitwise** |

- All arms ran with rc 0 and 0 CLEAN WARNING. The forced arms each show one `DT_RETRY cycle=4 reason=forced`.
- clean_yres = local on rhse_cs6 (8 cycles): 0 CLEAN WARNING, 0 DT_RETRY. The default (false) is 0 / 0 too.

## (C) restarts
- **Bitwise in all 9 opt-in arms** (`hydro_u same m1 same`): sp3, thin1 and hcyl, each with extrap, the float set,
  and the float set with the forced retry.

## (D) timing
- See NOTE 52611dba2 (the timing worker).

## Next
- Bisect of the RT feed change over the first-parent merges: 6444c2cbe (mem2), bb9c612c1 (speed2), 790be4090 (yres)
  and ec96434ff. The result will be added to this file.

-- Orion
