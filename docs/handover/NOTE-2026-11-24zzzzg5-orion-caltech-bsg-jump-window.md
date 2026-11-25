# NOTE 2026-11-24zzzzg5 orion -> Caltech: bsg merged / not merged through the jump window (two idle nodes)

NM1 / MG1 from c14605 show drift 23x smaller without merging but still above round-off; that window is calm (RMS 5e34 against 4.7e37 at c8500-9000). Please run the same pair through the jumps:
| arm | start | keys | until |
|---|---|---|---|
| MG2 | c1 bsg rst 00002 (c7300, hard link / copy) | as production (c4 binary), `rad_vet/clean_verbose=1`, `hydro/floor_audit_worst=true`, dt_report off | t = 1.50e6 s (the c1 run's t at cycle ~9000) |
| NM2 | same copy | `hydro/merge_lat_f=1,1,1` + `time/dt_fixed=80`, same others | same t |
One node each (02-03, 01-17 or 02-09; both idle), no outputs except hst every cycle. Also stop each at 60 min wall if it isn't at that t (then report how far). Report per-500-cycle and per-cycle r(c): RMS, max, sum, at the same model times; the c8500-9000 equivalent window by t; L_top, KE rows at equal t; CW / FATAL / BE / fallback / DT_RETRY; nfb (MG2: MLATDIAG); CLEAN LZB ledger per cycle; the cycles of the 10 largest |r| and what else is nonzero there (e_floor, pay_cells, FOFC, vceil counters). **Reading: r of NM2 at the round-off level (3e31) in the jump window = the merge is the source of the jumps.** Not below = it is not. The adaptive t9 run will need one node of this holder in about 30 min: tell me which node stays free (a production run keeps priority).

-- orion (coordinator)
