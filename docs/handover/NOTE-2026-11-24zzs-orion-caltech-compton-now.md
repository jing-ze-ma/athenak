# NOTE 2026-11-24zzs orion -> Caltech: USER: run the Compton arms NOW: pause bsg (not t9l62) to free a node

The USER wants the Compton test now, not after bsg's tlim. Pause **bsg** (it has ~17 h of slack against 24 h). Pausing t9l62, the critical star, would cost its own goal.
1. bsg: clean stop now (writes an rst; keep it). Its node 02-11 (7333073) becomes the test node.
2. Build caltech-t9hot-compton-f8e55fee7 (CUDA 12.9 + g++ 13, he_star_vet, separate dir; tiny gate without vet_gd_ls, pass 1 = 3.686e-11), if not built yet.
3. Arms (NOTE 24zzr), from an ERUPTION-phase t9l62 rst. The fresh T2 run isn't there yet, so use the rc3 run's `rst_keep/t9l62.final_1124.rst` (t ~5.7e5, inside the eruption; T2 restarts rc3 rsts with every key named + vet_gd_ls=0). Args: the T2 t9l62 production args (ia123rc4_decks/IA123_t9l62_restart.args, every key named) + `hydro/dt_report=10`:
   - C0: T2 binary (athena_ia123rc4v_262affd6f), key off, 300 cycles: the baseline on this rst;
   - C0b: compton binary, key off, 50 cycles: hst bitwise vs C0 (GPU key-off check);
   - C1: compton binary, `rad_vet/compton=true rad_vet/compton_report=50`, 300 cycles;
   - C2: C1 + `time/cfl_number=0.15`, to C1's end time (time-error yardstick);
   - C1r: C1 keys, 100 + rst + 100, bitwise vs C1.
   Run C0 + C1 at once if memory allows (2 GPUs each), then C2 / C0b / C1r.
4. Post per arm: dt every 10 cycles, DTREP / DTREG lines, `<rad_vet> compton` report lines, s/cyc, CW / FATAL, hst tail (L_top, L_in, E_tot, M, L_z).
5. Then resume bsg from its rst (same binary / keys) and post the relaunch line. Then CSHERO (24zzn) on that node only after bsg's tlim.

-- orion (coordinator)
