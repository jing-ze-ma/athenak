# NOTE 2026-11-23c orion -> Caltech: USER - prepare FRESH bsg, m35mms and t9l62 runs on the interim clean scheme IA9, cfl 0.3

USER (11-23 ~01:40 PST): prepare to simulate **bsg, m35mms and t9l6.2 from the start (fresh, NOT from a restart)** with the new
interim binary, using the ORIGINAL initial conditions (the plain bsg / m35mms / t9l62 decks), **not the rebuilt i1_* / ic_star ones**,
at **time/cfl_number = 0.3**.

## Binary
- GitHub jing-ze-ma/athenak, branch **`caltech-ia9-77551e823`** (Orion interim head IA9 = 77551e823).
- Build with **CUDA 12.9 + g++ 13** (not CUDA 12.2.1 + g++ 11.5: it miscompiles the ray closure on H200).
- Gate first: the tiny BSG smoke you used before; post its pass-1 numbers and the binary md5. (The old 4.020e-06 reference was for an
  older head; post what IA9 gives, we compare on our side.)

## What IA9 is (all of these are the code DEFAULTS on this branch; name them anyway in the decks)
- time scheme: rad_m1/time2_tableau=sdirk2, rad_m1/clean_stag_xdonor_freeze=4, rad_m1/time2_k1_stiff=0, rad_m1/time2_cfl_post=true
- faces: rad_m1/clean_w_fse=true, rad_m1/clean_face_hctau=true, rad_m1/clean_face_jse=true, rad_m1/clean_face_jse_cell=false,
  rad_m1/clean_bc_area=true, rad_m1/clean_stag_upwind_central=true, rad_m1/vet_seg_plateau=false
- restart/ring: rad_m1/clean_rot_ring_fill_rst=true, rad_m1/clean_rst_mshare=true (only acts on old-format restarts)
- keep each deck's own vet_rotate_avg / vet_rotate_avg_ema lines as they are; geometry / problem / opacity keys unchanged.
- time/cfl_number=0.3 on the command line for all three.

## Please
1. Build + smoke, post md5 + smoke numbers.
2. Prepare the three fresh decks (bsg, m35mms, t9l62; original ICs), with the keys above named explicitly, cfl 0.3.
3. Post the exact command lines, node counts and expected wall for each; launch when ready (same holders/allocations you use for
   production), and post the first ~200 cycles as usual: CLEAN WARNING / stage fallbacks / dt_retry, dt, s/cycle, L_top/L_in.
4. Report anything unexpected at startup (FATAL key checks, nghost: the clean flux solve needs mesh/nghost >= 3 where an x1 line is
   split across MeshBlocks).

Validation on our side (old-production restarts, own-window criterion): R136 0.102 % (300 s window, p99 0.129), RSG 0.357 % (p99 0.523),
HG L_top 2.69 % (p99 2.92); deck D 2nd order (rising to ~1.9 at the finest pair), restart/dt_retry bitwise, energy exact.

-- orion (coordinator)
