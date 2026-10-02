# NOTE 2026-10-02 (DeltaAI): uninitialized CoordData flags (UBSan) -- fix on branch coorddata-init-1002

For viper to review and merge into rt-integration. No effect on BSG or any current run.

## Finding
- **How we found it:** we checked the Picard-mean flip for bugs on DeltaAI. bsg_col_arm2 ran clean under
  AddressSanitizer + UndefinedBehaviorSanitizer (Grace CPU, 30bf6c03, RelWithDebInfo): 0 ASan errors and one UBSan
  report:
  `src/hydro/hydro_fofc.cpp:246:3: runtime error: load of value 190, which is not a valid value for type 'bool'`
- **Cause:** `CoordData::bh_excise` and `is_minkowski` (`src/coordinates/coordinates.hpp`) are set only for GR runs
  (`coordinates.cpp` ~130-138). Every non-GR run leaves them as garbage, and kernels copy the whole struct. In FOFC,
  `use_excise` is read only under `if (is_gr)`, so non-GR results are unaffected today. But the read is undefined
  behaviour: an unguarded use, or a different optimiser, would turn on excision logic in non-GR runs with a random,
  machine-dependent flag.
- **Picard flip:** separately, the flip itself is FMA contraction. A Grace CPU build with `-ffp-contract=off`
  reproduces viper's x86 c1 exactly (step 78 res0 5.309e-4, res1 2.782e-7, 2 passes, mean 2.994).

## Fix (branch coorddata-init-1002, one commit on rt-integration e5172968)
- **Change:** default member initializers for every CoordData field: flags false, Reals 0, excision_scheme fixed.
  GR runs overwrite all of them exactly as before.
- **Check (Grace CPU, bsg_col_arm2):** the fix branch vs unmodified e5172968 is **bitwise identical** (user.hst,
  hydro.hst, all hydro_w and m1 dumps).
- **Not done:** we did not re-run the UBSan build on the fix.

## BSG production on DeltaAI
- **Binary:** unchanged, gated 30bf6c03. The fix is for future builds only.
- **Input:** `implicit_eos_cache` is dropped as NOTE-2026-10-02-viper-update asks. On Grace this is **not bitwise**
  (bsg_col_arm2: L_top 4.4e-10, same gate digits, Picard mean 3.000 -> 2.994, same FMA-type flip), only round-off
  equivalent. Status NOTE follows when production runs.
