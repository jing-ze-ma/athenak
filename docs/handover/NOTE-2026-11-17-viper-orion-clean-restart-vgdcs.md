# NOTE 2026-11-17 viper -> Orion (USER): clean restart bug (stale vet_gd source below vlat_icut), fixed; production not affected

USER asked to tell you. It affects your staggered branches too if they read vgd_cs with clean FV rays; please check.

## Symptom
- R136 with the clean scheme was not restart-bitwise, even on the push-16 line, with or without any memory key.
- Restart at cycle 487124 + 3 cycles vs a 6-cycle continuous run: v 8.1e-5, F 4.9e-5, E 5.7e-8.
- Gate decks and your RT feed deck restarted bitwise, which is why the harnesses never saw it.

## Cause (found with a full-state fingerprint, then raw dumps)
- At the start of the first step after the restart, every carried array was identical between restarted and continuous: u0, opac, fluid,
  t2k1, cln_mn/mc, cln_fvi, csk1.
- The first difference was in the vet_gd ray-source array `vgd_cs` (ln chi, ln S) read by the first clean FV sweep:
  - VetGdSource copies vlat_cs into vgd_cs only for shells i >= is + vlat_icut, then halo-fills the lateral band. That is the range the gd
    sweep reads.
  - The clean FV sweep reads the WHOLE column.
  - On R136 vlat_icut = 294, so shells is .. is+293 kept the source written by the run's FIRST build, when vlat_icut was smaller.
- Raw dumps show the restarted vs continuous difference in exactly i = 3..296 after the first sweep, all k and j, both ranks. vlat_cs itself
  is identical in both runs.
- So clean runs depended on where they last (re)started: a stale, run-start-dependent source in the deep interior.

## Fix
fork clean-cc-rstfix-viper **5eb7adf3f** (off clean-cc 5785ee51d), 2 files, 6 lines:
- when clean FV rays are on, VetGdSource copies and halo-fills vgd_cs over every shell each build;
- VetGdCsSource's halo range gets the same change on the cubed sphere;
- the gd/production range is untouched.

## Evidence (MI300A)
- R136 fixed: restart vs continuous is data maxrel 0, with identical CLEAN lines.
- Straight-run change from the fix (6 cycles vs 5785ee51d): v 4.4e-4, F 4.6e-5, E 6.2e-8, rho 1.5e-9. It is a bug fix; the old value
  depended on the run's start.
- 21 push-16 gate decks: byte-identical (vlat_icut is 0 there). Your feed deck: straight vs restart bitwise, before and after the fix.
- **Production is not affected:**
  - the production gd sweep, moments and candidates only read shells >= vgd_scut;
  - measured: p1118 on R136, restart vs continuous, is maxrel 0 in all final dumps.

## New diagnostic
84fe1b09f, the next commit on that branch:
- `<rad_m1>/dbg_rst_fp=N` now fingerprints every cln_* array and scalar, vgd_cs, vlat_cs and fluid u0/w0.
- `dbg_rst_fp_cycle=C` adds fingerprints around each linear solve of cycle C.
- `dbg_rst_fp_dump=true` writes vgd_cs to files.
- Bitwise off and on on 21 decks.
- Use: run straight and restarted with the same keys and diff the RSTFP lines of the same cycle.

## Requests
1. Check stag2 (clean-cc-stag2-*) for any clean-path read of vgd_cs (or other vet_gd arrays) below vlat_icut / vgd_scut. If it is the same
   code, cherry-pick 5eb7adf3f.
2. Add a deck with vet_gd + vlat_icut > 0 (an sp deck with a deep interior, e.g. a he_star-like column) and a mid-run restart to reg/. None
   of the current reg decks exercise vlat_icut > 0.
3. This fix will be in push 18 together with speed round 3 (row-chunk FV sweep, bitwise; a safer clean_fv_lag = extrap fallback rule).
