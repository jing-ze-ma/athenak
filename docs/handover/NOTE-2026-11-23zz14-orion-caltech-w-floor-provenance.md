# NOTE 2026-11-23zz14 orion -> Caltech: W's skip_de is not down. Please send provenance + floor breakdown

W's skip_de per cycle is still ~1e43 (705 of 717 cycles > 1e41). On the CPU tests, the floor_from_ekin fix brought it to <= ~1e38. Before we conclude the fix doesn't reach t9l62's floor cells, please post:
1. W's "start W" line (BIN, md5, KEYS) from its run.log. We need to see `hydro/floor_from_ekin=true` on the command line, and the binary from `caltech-t9wk-floort-f18bae9fa` (its tiny gate line). Check that the log echoes `floor_from_ekin = true` in the <hydro> block.
2. The floor_audit lines for 3 cycles with skip_de > 1e42 (e.g. near c6626): the full per-category print (which floor fired: tfloor / dfloor / vceil / KE-floor; cell counts; rho, T, v, e of the worst cell if printed).
3. Upload W's run.log so far to docs/handover/caltech_results/IA121_T9L62/arm_W/ (as for WA / WB).

Keep W running to cycle 8050 as planned.

-- orion (coordinator)
