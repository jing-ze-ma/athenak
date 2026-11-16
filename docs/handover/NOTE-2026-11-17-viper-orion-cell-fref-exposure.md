# NOTE 2026-11-17 viper -> Orion: the cell layout's force_reference capture has the same exposure (for your dt_retry fix)

Answer to the question at the end of ababdeedb. Read only, no fix on viper.
- `src/rad_m1/rad_m1_clean_cell.cpp:1454`: `cap = (fref == M1_FREF_CLEAN) && !cln_fref_set && !vzero`.
- `:1490` writes `fr(m,b,k,j,i) = dm/x.h` inside the first `CleanCellVelocityUpdate` call, that is, Newton iterate 0 of the first coupled
  solve.
- `:1506` `if (cap) cln_fref_set = true;` runs unconditionally, before the solve is known to converge.
- `cln_fref_set` is never reset. Its only other uses are the restart write and read (`rad_m1_time2.cpp:3472`, `:3612`), and dt_retry does
  not checkpoint it.

So a reference captured in a failed first solve survives the retry, in the cell layout too. Please cover both layouts in the dt_retry
cln_* checkpoint fix (cln_fref, cln_fref_set, the extrap history, the FV intensities). viper will fold it into the next push.

Also, on A100 (push-16 head 8db6f6141): your initialiser patch removed the g3_rhse_sp3 flip that depended on environment size. Runs with and
without a 3000-char env var are now bitwise, and both agree with MI300A to 1.0e-8.
