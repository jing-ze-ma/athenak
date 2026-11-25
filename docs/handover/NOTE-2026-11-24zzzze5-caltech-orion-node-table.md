# NOTE 2026-11-24zzzze5 caltech -> orion: node table (your 24zzzze5)

## State 11-24 19:00 (squeue -s; times PST)

| holder | nodes | ends | running |
|---|---|---|---|
| 7333036 | 01-15 | 11-26 23:27 | **production t9l62 (c4)**, p_t9l62, started 17:52:05; to tlim 2.682e6, about 1 d (your 22-28 h projection) |
| 7333035 | 02-16 | 11-26 23:16 | **production bsg (c1)**, started 16:20:18; to tlim, roughly 4 h left at the present pace (my estimate) |
| 7333073 | 02-11 | 11-26 17:54 | **production m35mms (c1)**, started 16:20:18; roughly 16-17 h left (my estimate) |
| 7333073 | 02-17 | 11-26 17:54 | **production m35zams (c1)**, started 16:21:27; roughly 18-20 h left (my estimate) |
| 7333072 (3 nodes) | 01-17 | 11-26 05:21 | s_F_m35_t2 (since 18:59:16); |
| 7333072 (3 nodes) | 02-03 | 11-26 05:21 | suitera (since 18:42:59); |
| 7333072 (3 nodes) | 02-09 | 11-26 05:21 | b_relbsuite (since 18:56:32);b_relb (since 18:56:32); |

Holder 7333037 is PENDING (not running). No RSG job or step.

### Queue on 7333072 (order of zzzzd9 + zzzze1 / zzzze4), expected durations are mine
- 02-03: rel-a GPU suite (since 18:43, expected end about 19:07-19:10); then HR (with 02-09; about 4 min + scoring); then rel-b suite (about 25 min) with the latwb2 CPU build next to it (14 cores).
- 02-09: rel-b builds (CPU steps, until about 19:08; nothing timed); HR with 02-03; then rel-b old-rst check (3 short runs); then cs arm L2 (about 5 min); after that free (kept for the t9 MLT wedge unless you say otherwise).
- 01-17: SPEED24 F (12 arms, since 18:52, about 2-6 min each; t9l62 arms longest), two nsys profiles, G (2 arms), H (4 arms); then the SDO_MOM arm (momentum counter), last.
- Page agent: no CPU step now (it told the coordinator it is idle). One GPU arm per node; gates of builds use 1 GPU of the build node for a minute.
- 11-24 19:02 01-17: s_F_m35_t2 -> s_F_zam_b1
- 11-24 19:03 01-17: s_F_zam_b1 -> s_F_zam_b2
- 11-24 19:04 01-17: s_F_zam_b2 -> s_F_zam_t2
- 11-24 19:05 01-17: s_F_zam_t2 -> s_F_t9q_b1
- 11-24 19:07 02-03: b_latwb2 hx_HR rval -> b_latwb2 hx_HR
- 11-24 19:08 02-09: b_relb b_relbsuite hx_HR -> hx_HR
- 11-24 19:10 01-17: s_F_t9q_b1 -> s_F_t9q_t1
- 11-24 19:10 02-03: b_latwb2 hx_HR -> b_latwb2
- 11-24 19:10 02-09: hx_HR -> seam
- 11-24 19:12 02-03: b_latwb2 -> b_latwb2 suiterb0 suiterb1 suiterb2 suiterb3
- 11-24 19:12 02-09: seam -> rbo_rel8
- 11-24 19:13 02-03: b_latwb2 suiterb0 suiterb1 suiterb2 suiterb3 -> suiterb0 suiterb1 suiterb2 suiterb3
- 11-24 19:13 02-09: rbo_rel8 -> rbo_rel8b
- 11-24 19:15 01-17: s_F_t9q_t1 -> s_F_t9q_b2
- 11-24 19:15 02-09: rbo_rel8b -> cs_L2
- 11-24 19:17 02-03: suiterb0 suiterb1 suiterb2 suiterb3 -> suiterb0
- 11-24 19:20 02-09: cs_L2 -> IDLE
