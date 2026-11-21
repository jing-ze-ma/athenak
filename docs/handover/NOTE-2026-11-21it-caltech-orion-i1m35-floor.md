# NOTE 2026-11-21it Caltech -> orion: i1_m35mms was stopped by the USER (re 3fceb1ac2) -- rerun FRESH with tfloor 1000 K

Not a wrapper/binary end: the USER (in our session, 01:17 PST) stopped i1_m35mms (scancel of its step at cycle ~9650, t 9.10e4 s) and ordered a fresh rerun with the temperature floor at 1000 K. Reason: its deck (I1P_m35mms, ICREL) has tfloor = 3.3177582140e+12 = 24,800 K while T_eff ~ 29-33 kK and x1max ~ 1.1 R_ph, so 6-9 % of the cells just above the photosphere (20 % at the top) sat at the floor (the original m35mms with 1000 K has 6 floor cells). Please do NOT relink the old run.

New production run: **i1_m35mms_f1k** (prod_fresh_1121b/i1_m35mms_f1k, hpc-sm-01-17, step 6133592.601, started 01:22 PST), deck decks_prod1121/I1P_m35mms_f1k.athinput = I1P_m35mms with ONLY tfloor = 133780565496.86447 (1000 K, same code-unit factor as the m35mms deck); B_j3 athena_j3_a5257b7ee (md5 f71a22b6), prodg.sh, KEEPRST 3, END 17:25. Old run dir + rsts kept (prod_fresh_1121b/i1_m35mms). Please add i1_m35mms_f1k to the PROD monitor in place of i1_m35mms.

-- Caltech
