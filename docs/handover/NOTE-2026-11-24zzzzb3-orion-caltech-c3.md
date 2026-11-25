# NOTE 2026-11-24zzzzb3 orion -> Caltech: c3 (ray-hole cure, final form): build and TEST only

Branch `caltech-c3-ce559f800` = 95813cf85 + T9LATE's cure (+42/-14 lines, 2 files; no key). Reproducer converges in 3 passes; beams flux exact; one open test (idort_beams_mb with 4 blocks along x1: face-momentum difference 3.5e-7, only with several x1 blocks, which the stars do not have); Orion gates still running. **No production start yet (the USER decides).** Build it now (CUDA 12.9 + g++ 13, own dir, tiny gate, md5, .prov; `he_star_vet`), and when the HG arms (24zzzzb1) leave a node:
1. **t9l62 fresh, 1 cycle, then 100 cycles** (candidate_1124 deck, all keys): clean_init passes + residual, CLEAN WARNING count (must be 0), FATAL / inf / BE / clamped 0, s/cyc, dt.
2. **bsg, M35ZAMS, m35mms preflights on c3** (100 cycles each, fresh) plus, for each: L_top/L_in, KE, M, E at cycle 100 vs the c1 production run's first 100 cycles (same deck, 95813cf85): the size of the change. Post it.
The running production stays on 95813cf85. Production first.

-- orion (coordinator)
