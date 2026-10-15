# TASK Raven -> Caltech (10-15 night): re-gate the dt_retry fix (fin51-r-fix b16b66db) on the E1 + He giant seeds -- TEST ONLY

Raven worker DRB (/raven/ptmp/jinma/dtretry_bug_1015) found and fixed a dt_retry restore bug present in EVERY dt_retry binary
(5b7e1bc5, 21c4839d, 65c249b3, 1bd0c2ce, 5ab8c87e): at a vet_gd rotation-block start, the probe solve reads m1.esrc, which the failed
attempt left at w(t) (esrc not checkpointed; pgen static hs_mw_ not restored) -> non-bitwise after a retry, only while the MLT ramp is
active. Also: halo-cache records restored by key but not contents (multi-rank hazard), Linhao lbath_f not checkpointed (the L_in diff).
Fix branch **fin51-r-fix @ b16b66db** (= 5ab8c87e + 63919cdd + 0da548c5 + comment), pushed. Straight runs (no retry) bitwise to old builds.
Raven could NOT reproduce your post-ramp mid-block P!=B / C!=F (e1 gate 4552183): every Raven combination outside block-start-during-ramp
was bitwise. Please re-gate:
1. Build b16b66db nofma HOPPER90 (fin51 recipe, builds_fin51); record md5.
2. Re-run the rol2u-1015 gate (gate_cr2.sh, same seeds, debug QoS) with BIN=b16b66db for B, C, L1/L2, F (P stays the production binary):
   E1 (agcar3d.00052 seed) and He giant (hegiant.00060 seed). Report GATE-BASE / GATE-RESTART / GATE-RETRY.
3. Extra E1 arm F_all = F + `time/dt_retry_views=all` (checkpoints all Views): if F_all = C but F != C, the remaining cause is View
   state outside the minimal set; if both differ, it is non-View state. Report the first differing hst row for each.
No production change: agcE1r (5ab8c87e) is safe to continue (post-ramp: the esrc path cannot fire; any remaining diff is round-off grown
by chaos). If the re-gate passes, the coordinator will propose switching agcE1r/hgR2w to b16b66db at a link boundary (user decision).
