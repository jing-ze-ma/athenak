# NOTE 2026-11-21iw Caltech -> orion: production switch DONE (ki/kj/kl, USER-confirmed)

Binary B_fdef = athena_fdef_03b027502 (caltech-test-fdef 03b027502, CUDA 12.9 + g++ 13, md5 b7f21ec10fb735ecbe9dc325e9271d93), built in place on 02-11 (~5 min).
Smoke gate: pass 1 = 4.934e-06, NOT converged after 50, FATAL 0 (differs from 4.020e-06 only because fwb/fcone are fresh-run defaults now; not the trap).
Keys: rst keys + clean_conv_max 1e-4 + clean_face_opac harmonic + clean_stag_opac_update true + clean_stag_flim false + clean_stag_fwb true + clean_stag_fcone true + cfl (no damp, no gate).
Each run was stopped right after its latest complete periodic rst (B_j3 has no graceful stop), rst hard-linked into <run>/rst_backup_k9/:

| run | node | restart from | relaunched | cfl | first ~50 cycles |
|---|---|---|---|---|---|
| t9l62 | 02-09 | t9l62.00031 | 14:24 | 0.9 | CW 0, FATAL 0, dt 166 (was ~55-70 at cfl 0.3) |
| t9r | 02-03 | t9l62.00033 | 14:26 -> 14:37 | 0.9 -> **0.6** | at 0.9: 2 CW (c26418, c26420, stage fallbacks) + hydro dt COLLAPSE 134 -> 17 / 149 -> 24 => fallback rule: restarted at cfl 0.6 from the same rst at 14:37 |
| m35mms | 02-10 | m35mms.00015 | 14:26 | 0.9 | CW 0, FATAL 0, dt 17.7 (was ~7) |
| i1_m35mms_f1k | 01-17 | m35mms.00014 | 14:26 | 0.9 | CW 0, FATAL 0, dt 19.8 (was ~7.5) |

Wall END 17:25 (holders end 17:46). Logs: prod_fresh_1121b/<run>/log.N (newest). Note: t9l62 was first launched by mistake on 01-17 (env override) at 14:24 and killed within ~1 min during its restart read; relaunched on 02-09.

-- Caltech
