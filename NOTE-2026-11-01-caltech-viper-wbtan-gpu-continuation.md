# NOTE 2026-11-01 Caltech -> viper (USER): roche_wb_tan MHD continuation on GPU; g_t now stored in the restart file (wbtan-etafix-1101 @ 4029b621)

Reply to 2ca560ea / 794ba429 / 8a7277d7.

**GPU results (Caltech H200, 02-03), MHD r1m-type decks, restart pair 300 vs 150+150, compared from `<par_end>` on:**

| build | deck | result |
|---|---|---|
| 2a6b5e9f (integ, etafix + xe_init_pending) | gates_1101/integ_b_nokeys (defaults, no STS) | DIFF, 158M of 393M payload bytes |
| 952d586b = 2a6b5e9f + 267a32d0 (face-EMF restore; same content as your af982e409) | same | DIFF, 136M bytes |
| 952d586b, roche_wb_tan = false (control) | same | SAME |
| **4029b621** (fin51 d3f878be + etafix + EMF restore + g_t in rst) | r1m3 = trc_1101 r1m_trc + wb_tan on (production keys: outflow1, valfven_noinfall_fac 4, rkg_dt_safety 0.5, rkg_lean, cs_seam_flux shared) | **SAME** (MHD) |
| 4029b621 | w121r1_1101/hydro/r1h | SAME (hydro) |

- **E_B check (integ_b_sts deck, STS):** e1ecd210 spikes 1.797e31 -> 8.19e33; 2a6b5e9f stays at 1.797e31 through t = 32.
- **Your CPU pass on af982e409 vs our GPU diff on the same content:** GPU-only. Checksums of every array the probe touches (u0 w0 wder wtemp wbq0 uflx b0 bcc0 e3x1/e1x2/e2x3 eta_b ad_coef u1) are identical before and after the probe, at the fresh start and at the restart, and the restart g_t equals the fresh g_t (same checksum). compute-sanitizer initcheck: 0 errors (fresh C16 3 cycles; C32 restart 2 cycles). The 1-cycle restart diff is concentrated in the TOP cells (deep i < 10 identical). Root cause of the GPU-only difference with a probe at a restart: NOT found yet (state outside the checked list, e.g. a conduction/RT cache the probe's operator call refreshes, is the next suspect).
- **The fix we recommend (coordinator plan): 4029b621 on wbtan-etafix-1101** -- g_t is written to the restart file and read back; no probe at a restart:
  - generic per-cell pgen restart slabs (kPgenSlabRstMagic, behind the am_vector spin block; pgen.hpp pgen_rst_slabs / pgen_rst_slabs_in); unset = not one byte written, every old file readable;
  - wbtan_common.hpp: the probe points pgen_rst_slabs at g_t once valid; wbtan::LoadFromRestart on a restart (log "g_t read from the restart file (no re-probe)");
  - a file without g_t (old files, key switched on at a restart) re-probes with a log line: tested, restart of a d3f878be file runs and re-probes;
  - key off: MHD r1m C16 40 cycles bitwise vs d3f878be; hydro cs+roche state bitwise (user.hst Efloor columns: GPU atomics, as before).
- **Please:** cherry-pick 267a32d0 (if not yet) and 4029b621 into integ-test-1101 (the integ branch already has d0d1a738 as 1708c7a8 and your xe_init_pending commit; 4029b621 touches pgen.hpp/pgen.cpp/restart.cpp/outputs.hpp/wbtan_common.hpp/deep_hot_jupiter_rt.cpp). We keep hunting the GPU-only cause separately; it does not affect runs with 4029b621 (no probe at restart).
