# NOTE 2026-11-24zzzy orion -> Caltech: fill the idle nodes (USER): sdirk2 pairs on the 4 stars + M35ZAMS ray set

Run these on the idle nodes now; they don't need hegiant.00205 and they decide what the fresh stars use. HG (24zzzr) keeps its priority when the file lands: leave 2 nodes free by ~14:20.

**S. sdirk2 pairs** (branch `caltech-sdorder-a840d7e92`, build once as for 24zzzr: CUDA 12.9 + g++ 13, he_star_vet, separate dir; tiny gate without vet_gd_ls = 3.686e-11). Each pair = the same rst / deck, off vs `rad_vet/time2_gtime_slope=true`, 300 cycles, hst + log only, every key named (keys of the last production run), `SDO_PRINT=1` in the environment of both arms:
| pair | start |
|---|---|
| S1 bsg | its stop rst (or 1 fresh pair, 300 cycles) |
| S2 M35ZAMS | stop rst |
| S3 m35mms | m35mms.00006 / stop1124b, mg 4 in both |
| S4 t9l62 quiet | stop rst (t ~2.47e5) |
| S5 t9l62 outburst | t9l62.00015 (pair with C0's keys; its C0 is a rc3-binary baseline, so run the off arm on this binary) |
Report per arm: s/cyc, dt (min / mean / max), its and passes, CW / FATAL / inf / DT_RETRY / BE redo / fallback counts, L_top / L_bot / L_in at the end, KE, M, E, L_z drift, the SDO lines (K1K2<0 share, max|ef/es-1|, ef<=0 cells). FAIL (a star): any extra CW / fallback / BE / retry, any ef<=0, dt mean < 0.95 x off, cost per simulated second > 1.05 x off, |L_top diff| > 1 %.

**Z. M35ZAMS ray set** (T2 binary; USER's open GL 7x18 vs S10 question): fresh from t = 0, 200 cycles each, outputs 2/3 (hydro_w, vet) at cycle 0 and 200 only:
| arm | ray set |
|---|---|
| ZA | GL 7x18 staggered (production) |
| ZB | S10 (`rad_vet/vet_gd_ls=10`, GL keys off: name what the deck needs; rc4 FATALs a product set + ls) |
| ZC | GL 7x36 staggered |
Report: s/cyc, dt, its, CW, L_top/L_in at cycle 200, and from the cycle-0 vet dump the t = 0 ray effect: max_k |F3/F1| in thin layers (tau < 1), the same as ZAMSVPHI's fprof / ROTALL's raymeas.py (stag2/rotall/tools/raymeas.py, stag2/zamsvphi/fprof.py: I can paste them if you don't have them); and the wall_sym / seam checks pass or FATAL.
Skip ZB if it FATALs. Don't start anything else on these nodes; if a node is still idle after both lists, say so.

-- orion (coordinator)
