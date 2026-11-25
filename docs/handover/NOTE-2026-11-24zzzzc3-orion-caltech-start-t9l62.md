# NOTE 2026-11-24zzzzc3 orion -> Caltech: START t9l62 on c4 (USER go: "start t9 once the cure passes its gates")

Gates passed: t9l62 fresh 1 + 100 cycles on c4 (your de9dd1f16: init 3 passes, 1.235e-9, CW 0, FATAL / inf / DT_RETRY / clamped 0, one BE of the first step), Orion full suite 36/37 (the miss is the USER-accepted beams L1 re-pin 0.177), thinf = false byte-identical.
1. **Start now, production c4:** binary athena_c4_a000d803c (md5 3ae130e9, checked by the launcher before srun), deck candidate_1124 t9l62 (md5 f2e3b57d) copied unchanged to a NEW dir `prod_c4_1124/t9l62`, fresh, no command-line keys beyond what the deck names, deck outputs and rst cadence unchanged, chain restarts only as for the other stars. Deck diffs: say none.
2. **Node:** one node, 4 H200, 4 ranks, in a holder with >= 34 h left (the projection is 22-28 h). Don't displace bsg, m35mms or M35ZAMS. If no holder has that much, say so before starting; don't start a run that must be cut. Tests on the chosen node: let a running one finish, then give the node up; the queued arms (zzzzb9 done, zzzzc2, zzzzb6, zzzzb5) use the other idle nodes.
3. Watch as c1 (watch.sh): status at c100, 500, 1000, then every 10 min; CW, fallback, clamped, FATAL, inf, DT_RETRY. FATAL or inf: clean STOP, no relaunch, post it. First rst hard-linked in rst_keep.
4. Post the start line (chunks.txt: node, holder, end, md5s), clean_init line, and the first 100 cycles (dt, s/cyc, L_top/L_in, floor_pay). 
5. The other three stars on 95813cf85 continue untouched. **Do not stop or restart them;** the bsg / M35ZAMS / m35mms c4 preflights (24zzzzb8 item 3) are still wanted and decide that later.

-- orion (coordinator)
