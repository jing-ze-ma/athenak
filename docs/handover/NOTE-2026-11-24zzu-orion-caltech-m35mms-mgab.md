# NOTE 2026-11-24zzu orion -> Caltech: USER OK: pause m35mms for a 2-vs-4 multigrid-level A/B (~10 min), then resume

m35mms's s/cyc rose 0.65-0.72 -> 1.17 (c5950, mg 2). The USER agreed to stop it briefly to test.
1. Clean stop of m35mms now; keep the final rst (hard-link rst_keep/m35mms.pre-mgab.rst). Same node (02-17), T2 binary (md5 7d7c4a71).
2. From that rst, production keys (every key named), outputs off (hst + log only, no rst writes), 200 cycles each, one after the other, with nothing else on the node:
   - A: rad_vet/clean_stag_mg_levels=2 (as production)
   - B: rad_vet/clean_stag_mg_levels=4
   - A2: repeat of A (noise check)
3. Report per arm, over cycles 50-200: s/cyc (mean and min), Krylov its and solves per stage, passes per cycle, dt, CW / FATAL / DT_RETRY; hst end values L_top, E_tot, Lz_tot of A vs B (differences are solver-tolerance level).
   Also say whether the 1.17 s/cyc itself is reproduced by A: if A shows ~0.7, the slowdown was from something else (neighbour load, outputs), and please say what (nvidia-smi clocks / other steps on the node).
4. Resume production at once from the pre-stop rst with mg 2 (the USER's uniform choice) and post the relaunch line. Don't switch to 4; the USER decides after seeing A / B.
5. Same question, cheap: what is bsg's s/cyc at its latest 300 cycles vs its first (0.69)? If it also drifts up, say so.

-- orion (coordinator)
