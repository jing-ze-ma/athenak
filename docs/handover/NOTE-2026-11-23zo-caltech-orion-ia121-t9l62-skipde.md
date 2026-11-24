# NOTE 2026-11-23zo Caltech -> orion: t9l62 IA12.1 after the transient + IA11 comparison (re 23zn "also please post")

caltech_results/IA121_T9L62/:
- prod_floor_audit_cw_6621_now.txt.gz: floor_audit, CLEAN WARNING, DT_RETRY and pos_floor lines for cycles 6621-6948 (~18:55).
- ia11_clean_t2s2_6450_6620.txt.gz: IA11 (prod_ia11_1123/t9l62/log.1) "CLEAN cycle N t2s 2" lines for 6450-6620. **IA11's log has no floor_audit lines**: floor_audit was off in that production run, so there is no IA11 skip_de to compare.

Reading of the production IA12.1 run:
- skip_de in [update] is CUMULATIVE: 7.548e42 at 6621, 7.5919e42 at 6681, 7.59258e42 at 6741, then flat (7.592582e42 from ~6800 to 6948).
  So the per-cycle increment has been < 1e40 on every cycle since 6621 and ~0 since ~6800. The transient has settled.
- Since 6620 there is 1 CLEAN WARNING (6637, t2s 2 stage fallback 8) and no new implicit_pos_floor.
- At 6900: t 4.7707e5, dt 9.29.
- The hourly skip_de-increment check (> 1e41 followed by a fallback -> post) is added to our watch.

Arms S / L on N2 node B (hpc-sm-02-17): queued. Build caltech-t9arm-932dbd534 -> tiny gate -> S and L (2+2 GPUs if memory allows, else S first). Results follow.

-- Caltech
