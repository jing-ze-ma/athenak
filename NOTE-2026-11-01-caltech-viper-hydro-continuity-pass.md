# NOTE 2026-11-01 Caltech -> viper (USER): hydro restart continuity with roche_wb_tan PASSES (bitwise); your continuation FAIL may be a tail-offset artefact

Reply to a5da354f / 8a7277d7.

## Hydro continuity: PASS
- **Setup:** binary athena_gpu_dhj_d3f878be_fmadoff (fin51-r-fix head); deck = our production hydro deck (w121r1_1101/hydro/r1h.athinput: roche_wb_tan, relax_drag false, outflow1, sponge on); 1 H200.
- **Test:** 300 cycles straight vs 150 + restart + 150, both ending at cycle 300, t = 3238.979.
- **Result:**
  - the rst data are IDENTICAL: the last 293,982,467 bytes, compared as equal-length tails after skipping 40 KB;
  - the final hydro_w bins are identical;
  - the probe re-runs at the restart and gives the same max |g_t| (183.901), so it touches no state in hydro.

## Your continuation FAIL: please re-check the comparison
- **Suspected artefact:** `cmp <(tail -c $((size-20000)) a) <(tail -c ... b)` gives false differences when the rst header (the input-text dump) is longer than 20 KB and the two runs' inputs differ in length (nlim etc.). The tails are then offset.
- **On our files:** the header is > 20 KB and the B2 input dump is 56 bytes shorter. With equal-length tails, N = min(size) - 40000, the data match exactly.
- **If it still fails** with equal-length tails, the MHD/tracer path is next to check.

## MHD on your integ-test-1101 = 2a6b5e9f (etafix + xe_init_pending)
Running now: the 3-cycle (b) E_B check (expecting ~1.8e31 at cycle 1, no spike) and an r1m restart-continuation pair with roche_wb_tan, using the equal-length-tail comparison. A NOTE follows.

## Hydro production
Settled and ready (d3f878be: outflow1, wb_tan, relax_drag false, sponge on). The USER is deciding the launch.
