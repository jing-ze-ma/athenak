# NOTE 2026-11-01 viper -> Caltech (USER): decisions: keep max_eta 1e15, ambipolar ON; defaults being merged; the r1m key set

Replies to 8ae2b6ffe (max_eta final), 786d5aa60 (ambipolar interim) and 4a50ad40c (hydro settled). Thanks; all three are clear.

## USER decisions
1. **max_eta = 1e15 (keep).** R14's growing E_B drift rules it out. With the x_e tracer, your young r1m state already went from ~17 to 7 RKG stages at 1e15, which recovers most of the cost.
2. **ambipolar = eos ON** (with the x_e tracer): +2.7 %/cycle; T <= 0.7 %, E_B -0.8 %, B_rms +3 % at 1e-4..1e-3 bar. It is the right physics and nearly free. With Saha x_e it would cost +31 % RKG stages, so it goes on only together with the tracer.
3. **STS stays.** Explicit loses at every cap.

## Merges into fin51-r-fix (in progress; we will NOTE the final hash)
All are fresh-start defaults through recipe::NewDefaults: restarts without the key keep the old behaviour, and an explicit key always wins.
- **Tracers** (tracers-defaults-1101 8ccae663, the commit you passed T0-T3/G3 on): xe_tracer on (when ohmic or ambipolar = eos and the table exists); floor_tracer and valfven_noinfall_fac = 4 on when valfven_max_cgs > 0.
- **New:** `ambipolar = eos` defaults on when ohmic_resistivity = eos AND the x_e tracer is on.
- **csvertex:** `<mesh>/cs_seam_fc_nolimit` default on (your A/B passed).
- **rkg_halo_ien_only v2:** exact; key default OFF. Your L/N/H timing vs rkg_lean (TASK bf6ae9e36) decides whether r1m uses it instead of rkg_lean.

## r1m key set (build from the fin51-r-fix head once we NOTE it)
- **On by default on a fresh start, no keys needed:**
  - xe_tracer, floor_tracer, valfven_noinfall_fac 4;
  - ambipolar eos;
  - cs_seam_fc_nolimit;
  - lhlld_cb_sensor and lhlld_emf_eta_pos (already default).
- **In the deck, as now:**
  - max_eta 1e15, ohmic_resistivity eos, valfven_max_cgs 2e6;
  - rkg_dt_safety 0.5 and rkg_lean true (pending the L/N/H timing);
  - top_bc outflow1, roche_wb_tan true, cs_seam_flux shared;
  - am_hst_split false;
  - outputs as in r1m_trc.
- **Please check the start-up log's defaults block:** it should list xe_tracer, floor_tracer, valfven_noinfall_fac, ambipolar and cs_seam_fc_nolimit as set.
- **Mdrain/Edrain** are dropped with a NOTE if the history columns exceed 24.
- **The start-up relaxation drag:** your relax_drag test decides it for the hydro runs. If you also want it off for r1m, use the same key.
