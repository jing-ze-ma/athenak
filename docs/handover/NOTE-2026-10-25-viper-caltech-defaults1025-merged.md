# NOTE viper -> caltech (reply to TASK bc221a716, step 2): defaults-1025 gated, fin51-r-fix = 83d97f8c9

**fin51-r-fix fast-forwarded c60714bd1 -> 83d97f8c9** (defaults-1025).
- **How it was gated:** on viper, CPU, login node. No GPU confirmation: apudev is down and apu is full.
- **Details:** viper `/viper/ptmp2/jinma/defaults1025/GATES.md`.
- **Decks:**
  - the HG production deck shrunk to 704x16x16 (4 ranks, 40 cycles);
  - the gateA grey sp wedge (300 cycles), and the same wedge with `<mhd>`;
  - the Cartesian 2-D He slab with hesdirk2.
- **Comparison:** history rows + bin data arrays, bitwise. Whole-file md5s cannot be used, because the headers carry the input text.

1. **PASS.** New binary, no keys == old binary with the production keys named, on HG and on the wedge.
   New vs old with no keys differs, so the changed path is exercised.
2. **PASS.** New binary with the old values named == old with no keys, on HG and on the wedge.
3. **PASS, with notes.**
   - **MHD:** posfix stays 0 and there is no FATAL. The sp MHD deck takes the new work_fix + extrap2 defaults, as designed
     (bitwise == old with those two named).
   - **Cartesian:** CHANGED by the work_fix default, as designed (new == old with only `implicit_work_fix=true`, bitwise).
     Every rad_m1 deck, Cartesian included, now gets the work fix unless it names `implicit_work_fix=false`.
4. **PASS.** A restart with an input lacking the keys: vet_col = lag, and work_fix/posfix take the new defaults (bitwise == old
   with them named).

**Nuances:**
- **vet_col on restart:** every sph run writes the vet_col value into its rst (`time2_vet`), so a restart of a current sph run
  keeps that run's value, e.g. predict.
  - It falls back to lag only for rst files from binaries older than the key.
  - Same behaviour as before.
- **posfix:** the re-solve path never triggered on these decks (0 re-solves).
  - The default of 4 is shown to be selected (bitwise == 4 named).
  - The re-solve itself is covered only by your production runs.
