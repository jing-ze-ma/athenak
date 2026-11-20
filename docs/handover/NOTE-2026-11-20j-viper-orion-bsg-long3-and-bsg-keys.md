# NOTE 2026-11-20j viper -> orion (USER): BSG long3 (vstar off + full tensor) result, and the fresh BSG run uses Linhao's GM, μ and floor

## 1. Fresh v20.1 BSG run (update to NOTE-2026-11-20i): use HIS values (USER 11-20)

- **μ = 0.62** (Ma+2026 text / his deck); `temp_unit_kelvin` = 7.4568864364e-9.
- **GM = 2.62502e27** (19.78 Msun, his deck value 260.795 code).
- **dfloor = 6.15245e-14** (his athinput.bsg_star.rusty: 1e-9 x density unit). Keep his floor. It is far above the
  1e-16 of the T9/T19 corona tests, so this BSG has no near-empty cells.

The staged `fresh_prod_1120/bsg/bsg.athinput` now has μ 0.62 again; NOTE-20i had staged 0.614. README and MD5SUMS are
updated. For information only: his profile alone balances at μ 0.614 with its own M(r), but the USER keeps his values
so the setup matches his run.

## 2. BSG long3 result (TASK-2026-11-07-bsg-novstar-raven; this report was never sent)

The run is `prod_lhxo` with two changes: radiation Doppler/work velocity = fluid velocity instead of vstar, and the full
tensor on. Everything else is the same, including he_bc_linhao_fis = iter.
- Fresh start to 57.4 d on Raven. Binary athena_he_a100_nofma_c_698b8241_fmadoff_clean (prod-1107v = prod-1107 +
  vstarfix-1106), not 3b1d3dc19 itself.
- Two links (32025175 to 52.2 d, then 32037236), ended on time limit 11-09; no FATAL. STEP-0 full-tensor gate passed.
- Run dir /raven/ptmp/jinma/bsgnv_1107/long3. Dashboard page https://claude.ai/artifact/2Mv9Qe8eNmAw94kCYfJeHw.

From the two `bsg3d.user.hst` files (L_top relative to its mean in each window):

| | long3 (vstar off, full tensor) | prod_lhxo (vstar on) |
|---|---|---|
| L_top min-max / mean, 30-57.4 d | 0.968-1.057 | 0.757-1.193 |
| L_top rms, 30-57.4 d | 2.2 % | 9.3 % |
| L_top rms, 10-30 d | 7.4 % | 11.6 % |
| KE_int mean, 10-57.4 d | 1.2e42 erg | 1.3e43 erg |

(Ma+2026: L range 0.91-1.09.)

**Verdict:** turning vstar off does not raise the variability toward Ma+2026. long3 has about 10x less KE_int and a
quieter light curve than prod_lhxo. The change also folded in the full tensor, and no tensor-off arm was run, so the
two effects are not separated.

**Caveats:**
- The window-mean normalisation includes the slow decline of L_top, which is larger in prod_lhxo.
- These prod_lhxo amplitudes do not match the "~0.98-1.02" quoted in the TASK. Please check against your own
  prod_lhxo analysis.

-- viper (for the USER)
