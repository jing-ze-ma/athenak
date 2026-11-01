# NOTE 2026-11-01 viper -> Caltech (USER): more of your production keys become fresh-start defaults; rkg_lean stays (thanks for the L/N/H timing)

Replies to 3ff4c8dd9 (L/N/H), f106d031e (vet_gd check) and d07e73225.

## Your L/N/H result
- Accepted: rkg_lean stays for r1m.
- rkg_halo_ien_only v2 is still being merged with key default OFF. It is exact, but it only pays on multi-rank non-lean runs, and valfven_noinfall_fac forces the full slot set anyway.
- The vet_gd tip finding is noted. Good that HG and RSG are clean.

## Our final x_e tracer run (s10 developed state, 2 rot)
- Clean, and the night F is stable (1e-6..1e-2 bar: 3900 / 4200 / 620 / 60).
- RKG stages fell from 34 to 22 per cycle by the end, so the tracer saving appears on developed states too, only more slowly than on your young r1m. Our earlier 3-6 % figure was an early snapshot.

## USER: additional fresh-start defaults (worker running; NOTE with the hash when merged)
Same rule as the others: recipe::NewDefaults on fresh starts only. A restart without the key keeps the old behaviour; an explicit key always wins.
1. **STS:** use_rkg_sts = true, rkg_lean = true, rkg_dt_safety = 0.5, whenever a resistive diffusivity is on. rkg_lean only where valid (uniform single-level cs).
2. **dhj top and tide:** problem/top_bc = outflow1, and problem/roche_wb_tan = true when problem/roche is on.
3. **Seams and history:** cs_seam_flux = shared, and problem/am_hst_split = false.

NOT a default: max_eta. Keep 1e15 in your decks.

## Together with the merges in flight
- xe_tracer, floor_tracer, valfven_noinfall_fac 4, ambipolar eos (with the tracer), cs_seam_fc_nolimit.
- A fresh r1m-type deck will then need almost no special keys.
- Your existing decks stay valid: explicit keys win.
- Please check the start-up defaults block on the first launch.
