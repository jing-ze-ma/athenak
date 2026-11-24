# NOTE 2026-11-24zzh orion -> Caltech: GO: build T2 (the ONE binary) and start the 4 stars FRESH on it

**Binary:** branch `caltech-ia123rc4v-262affd6f` (pushed). It's the release trunk + rad_vet trim + rc4 + DFLFIX + floorcap_lz_spin + rad_signal_speed_phase + dt_floor_courant + ROTWB (inert without rotation).
- Orion gates: key-off suite 22/22, 77 dirs bitwise; every production key group and pair restart-bitwise; never-read 10/10 empty; R136 / RSG hand-overs rc 0 / 0 CW; tiny gate PASS 3.686e-11 (md5 on Orion e56c939c).
- Build: CUDA 12.9 + g++ 13, he_star_vet, separate build dir. Tiny gate deck WITHOUT vet_gd_ls (rc4 FATALs a product set + ls). Expected pass 1 = 3.686e-11.

**Decks:** docs/handover/ia123rc4_decks/ (tile/ and notile/; MD5SUMS). They name mg_levels 2, dt_floor_courant 0.6, rad_signal_speed_phase, fofc + fofc_retest, floorcap_lz + _spin + mesh_vol, seam keys, clos1_budget + bud_maxnorm. **Fixed GL 7x18 rays, no rotation:** ROTALL measured +12 % s/cyc for rotation, above the USER's cost bar.

**Order:**
1. Let the arms already running finish (300 cycles each). **Drop the ZAMSFAST arms** (24zzg): T2 has spin + dt_floor_courant, so the fresh M35ZAMS run covers them. Don't start new T9HOT2 / T9FAST arms. T9HOT P1 / P2: finish if started, else drop. I'll ask for one cfl 0.15 yardstick on T2 later if needed.
2. **Start fresh, each on its own node, each star's best layout** (as the rc3 chunks: bsg / m35mms / t9l62 tile decks; M35ZAMS the Z3 layout):
   - **t9l62 → 7333035 or 7333036** (longest holders; ~52 h left);
   - **M35ZAMS → the other of 7333035 / 7333036**;
   - **bsg + m35mms → 7333073's two nodes** (bsg ~5 h; m35mms ~20-23 h, inside 11-26 17:54).
   - The first 300 cycles of each fresh run are the gate: post a 300-cycle watch per star (s/cyc, dt, binder from DTREP if named, Krylov its, CW / FATAL / inf / DT_RETRY). Stop rule: CW / FATAL / inf / DT_RETRY only, never dt alone.
   - Name `hydro/dt_report=10` for the first 300 cycles only (the deck can't take it on a fresh command line unless it's named there, so add it to your deck copy), then relaunch without it at the first rst. Or skip it if that's awkward.
3. Then the usual 30-min status line + ETA table every ~2 h.
4. When bsg finishes (~5 h), its node becomes the test node (R136 on T2 later; ask me first).

-- orion (coordinator)
