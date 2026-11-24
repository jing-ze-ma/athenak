# NOTE 2026-11-23zl Caltech -> orion: t9l62 IA12.1 lines for cycles 6450-6620 (re 23zk)

Posted in caltech_results/IA121_T9L62/:
1. clean_warn_dtretry_6450_6620.txt.gz: all CLEAN / WARNING / DT_RETRY lines with cycle 6450-6620, plus the untagged DT_RETRY / floor_audit lines (3001 lines).
2. wbfdiag_floor_audit_6450_6620.txt.gz: the WBFDIAG and floor_audit subset (546 lines).
3. t9l62.user.hst.gz and t9l62.hydro.hst.gz: every row from the 17:47 restart (new dir prod_ia121_1123/t9l62) up to 18:18.

t9l62 keeps running (USER also said keep running and watch). Caltech now watches every 5 min with your 23zk criteria: dt < 0.01, NaN/inf/non-finite, FATAL, or "implicit_pos_floor created". On any of them we stop cleanly and post (USER approved stopping on these). N2 stays idle until your next NOTE.

-- Caltech
