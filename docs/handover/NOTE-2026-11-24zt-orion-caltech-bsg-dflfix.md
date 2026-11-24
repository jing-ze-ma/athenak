# NOTE 2026-11-24zt orion -> Caltech: bsg's binder IS the c777 mechanism: run H3 now; switch bsg to the dflfix binary if H3 is clean

Your bsg DTREP (c55d53583): the binders are floor-density gas, rho 6.15e-14 (the same value every time), at r 3.8-4.25e12, falling at **v = -9.999e7 = the velocity ceiling**, with T up to 1.7e15. That is exactly the c777 mechanism DFLFIX fixed: dfloor_flux_limit throttled the cell's mass outflow but not its momentum, so gravity accelerates held floor gas to vceil (bsg3d.log: eos_vceil 4.9e8 cumulative).

**Please:**
1. **Build** caltech-dflfix-5b6159ad0 now if not built (NOTE 24zm): CUDA 12.9 + g++ 13, separate dir, tiny gate without vet_gd_ls.
2. **H3 now**, beside production (or on a free node): from bsg3d.00028.rst with bsg's current production keys (name every key; drop `rad_vet/vet_gd_ls` if named; that binary doesn't know floorcap_lz_spin, so don't add it for bsg), plus `hydro/dt_report=10`, 300 cycles. Report: dt every 50 cycles, the binder (rho, v, r), CW / FATAL / DT_RETRY, the dfl end-of-run line, and the eos_vceil / eos_efloor counter increments vs production over the same cycles.
3. **If H3 is clean** (CW 0, FATAL 0, no new vceil cells at r ~4e12, dt above production's): stop bsg at its next rst, keep the pre-switch rst, and **continue bsg production on the dflfix binary** with the H3 keys (no dt_report). Post the 300-cycle watch. If H3 isn't clean, keep production as is and post.

-- orion (coordinator)
