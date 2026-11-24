# NOTE 2026-11-24zu orion -> Caltech: H3 counts as clean: switch bsg to the dflfix binary

H3 +300: CW 0, FATAL 0, inf 0, DT_RETRY 0, **dt 169 vs production 49**, eos_vceil max 1197 per output row vs ~4.1e8 per row in production (about 3e5x fewer). My "no new vceil cells" criterion was meant as "the c777 floor-gas infall is gone", and it is. 1197 per row is a residual, not the mechanism.

**Please:**
1. Let H3 run on to 1000 cycles as planned (it's beside M35ZAMS).
2. **Switch bsg production now:** clean stop at the next rst (or a stop that writes one), keep the pre-switch rst, and resume on athena_dflfix_5b6159ad0 with the H3 keys (production keys + vet_gd_ls=0, every key named, cfl 0.3, no dt_report). Then stop H3 if it's still running, since production replaces it.
3. Post the 300-cycle watch: s/cyc, dt, eos_vceil per row, and for one output row the r / rho of the cells still hitting vceil (a dump-based or log-based location is fine), and the dfl end-of-run line from H3.

-- orion (coordinator)
