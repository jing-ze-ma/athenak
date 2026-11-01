# TASK 2026-11-01 viper -> Caltech (USER): time the exact RKG halo saving (rkg_halo_ien_only v2) against rkg_lean on the r1m deck

The USER asks for this comparison. Short: cost runs only.

## What v2 is (fork/rkghalo-test-1101 @ 3503c1661; being merged into fin51-r-fix, key default OFF)
- `<mhd>/rkg_halo_ien_only = true`: RKG super-stages >= 2 exchange only u0(IEN) plus the slots the previous stage's ConToPrim changed. On W121 that is {IDN, IEN} every stage, because the density floor fires; never all slots.
- The changed set is found by a per-slot device reduction, OR-combined across ranks. It is EXACT: bitwise identical to key off.
- **Viper MI300A, s10 non-lean RKG:**
  - GPU bitwise PASS on s10, b100 and s10 + ambipolar;
  - cost 0.2395 -> 0.1989 s/cycle (-17 %).
- It is FATAL together with rkg_lean (they replace the same exchange).

## Question
r1m uses rkg_lean, which skips the intermediate C2P, the U halo exchange and dt, but is NOT bitwise. Is lean still meaningfully faster than the exact path (non-lean + rkg_halo_ien_only), given that lean differs at the cfl-twin noise level?

## Arms
- Binary: 3503c1661 (or the fin51-r-fix head once we NOTE the merge).
- Start: one r1m rst (e.g. off_0p50.rst). Production keys, rkg_dt_safety 0.5. Same max_eta as your current best (R15, plus R14 if cheap).
- Run interleaved on ONE dedicated H200, 3 repeats x ~300 cycles each:
  - **L:** rkg_lean = true (production);
  - **N:** rkg_lean = false;
  - **H:** rkg_lean = false + `mhd/rkg_halo_ien_only=true`.
- Report:
  - s/cycle per arm and RKG stages/cycle;
  - the `rkg_halo_ien_only:` histogram line from H.
- Optional: L vs H difference after ~0.5 rot (winds, T, B_rms), to see whether lean's non-bitwise noise stays at the cfl-twin level.

## Decision (USER)
- If H costs within ~10 % of L, we would recommend the exact path (H) for r1m.
- Otherwise lean stays.
