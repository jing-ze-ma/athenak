# NOTE 2026-11-21jb Caltech -> orion: B_p1 production switch DONE (ky/kw/lb; USER confirmed in our session 16:33)

B_p1 = athena_p1_bf4f917d4 (caltech-prod-bf4f917d4, CUDA 12.9 + g++ 13, in place, ~4 min), md5 d38d411c5f29f70370efc2061f060b1b.
Smoke gate: pass 1 4.934e-06, NOT converged after 50, FATAL 0 (identical to B_fdef/B_fwbs class; fwb/fcone fresh defaults; not the trap).
Keys: rst keys + conv_max 1e-4 + face_opac harmonic + stag_opac_update true + stag_flim false + stag_fwb true + stag_fcone true + stag_fwb_stage true + stag_thinf_held false + clean_gas_newton_switch 4 + clean_gas_newton_switch_min 6 + cfl; env VGD_XFE=0. No damp, no xpat.

| run | node | restart from | cfl | extra | relaunched |
|---|---|---|---|---|---|
| m35mms | 02-10 | m35mms.00017 (last good, before the 16:16 FATAL) | 0.9 | | 16:37 |
| t9r | 02-03 | t9l62.00034 | 0.6 | hydro/tfloor 1000 K | 16:37 |
| t9l62 | 02-09 | t9l62.00035 | 0.9 | | 16:37 |
| i1_m35mms_f1k | 01-17 | m35mms.00018 (waited for it to close) | 0.9 | | 16:39 |

Each rst hard-linked into <run>/rst_backup_p1/. Wall END 17:25 (holders end 17:46). First cycles follow.

-- Caltech
