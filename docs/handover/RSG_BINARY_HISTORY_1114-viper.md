# RSG_BINARY_HISTORY_1114: which binary has been most stable for RSG E production (read-only audit, 11-14)

Sources (all read in place, nothing launched or edited):
- viper: /viper/ptmp2/jinma/prod_1107/rsg (run.cfg.<job>, cmd.<job>.txt, holder.log, run.*.log)
- Raven: /raven/ptmp/jinma/raven_gap_1105/rsgE (same layout, plus caltech_hst/ = Caltech prodH history to 6326.6 d)
- Caltech prodH: /raven/ptmp/jinma/caltech_1107/rsg_E_1020/{runH.cfg,prodH/inalloc.out} (START lines only; the Caltech run logs are not on viper/Raven)
- Counts are `grep -a` over each run log (fofc_report lines skipped); cycles, sim time and wall from the `elapsed= cycle= time= dt=` lines
  (script /tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/42a91416-7b38-4094-88eb-372bd1345dea/scratchpad/seg.awk).
  NA = "NOT ADMISSIBLE ... redone with backward Euler" (= BE redo; dt_retry lines track it 1:1). BE share = BE steps / all steps from the
  end-of-run hesdirk2 summary (only printed on a clean exit). Picard/cyc = mean of rsg.user.hst column 22 (the column never reaches 200, so it is not an NC detector).
- All times CET (machine clocks). Hardware: Raven = 24 x A100 (6 nodes), viper = 8 x MI300A (4 nodes); s/cycle is only comparable within a machine.

## 1. Segments (one row per production link)

| # | machine / job | binary (run.cfg snapshot / cmd file) | key change vs previous | sim t [d] start -> end (span) | cycles | wall h | s/cyc | mean dt [s] | NC | NOT POS | NA (BE redo) | FATAL | BE share | end |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 | Caltech 5691175 et al. | e04fbd50 vtxfix2 -> d85bfbce autobe2x -> 60f7ccda/e1bbc562/651892bc/12b7deef/3e0bd80a -> 631d1672 -> prod1102_6c76e282 | vtxfix restart from 1675 d (after the prodG 1695-1704 d vertex ejection; prodG also had a 562 d corner eruption) | 1675.0 -> 4505.0 (2830) | 258,855 (hst rows) | n/a | n/a | ~945 | n/a | n/a | n/a | n/a | n/a | binary switches only (inalloc.out) |
| C2 | Caltech | prod_23a55658_fmadoff (prod-1102 + 11-03 seam keys, `cs_seam_rho_guard=0` g0_1103) | seam_blockend=extend, merge_lat=exchange, vertex_3cell; guard 0 | 4505.0 -> 6326.6 (1821.6) | 167,443 (hst rows) | n/a | n/a | 940 | n/a | n/a | n/a | n/a | n/a | handed to Raven 11-07; Picard/cyc 8.5 |
| R1 | Raven 32026499 | athena_he_a100_nofma_c_698b8241_fmadoff_clean (prod-1107v) | r8 keys: **guard 0 -> 4**, cs_cross_mono, implicit_fb_local(4), opac_lag_passes 3, predictor_tcap 1.1, k1_redist 0.1; vet_gd_tangential (full tensor); layout_invariant | 6325.1 -> 6449.0 (123.9) | 10,750 | 4.72 | 1.58 | 996 | 0 | 1691 | 478 | 0 | n/a | srun rc=143 (external SIGTERM, not a code abort); holder held, FIXED -> R2 |
| R2 | Raven 32026500 | athena_he_a100_nofma_o_864d120d_fmadoff_p1109 (prod-1109) | + hydro/wall_wb_onesided, hydro/split_c2p_once (db4af9ce/p1108 never ran a link) | 6445.5 -> 6987.5 (542.0) | 48,000 | 22.50 | 1.69 | 976 | 1 | 9425 | 3127 | 0 | 13.2 % | rc=0 wall limit |
| V1 | viper 12796542 | athena_he_gpu72_p_4e9f602c_p12r (rstabe) | + **geo2 + EMA128 + 3-D rotation** (vet_gd_geo=2, gl_nmu=0, rotate_dt=8000, align, axis=3d, ema step), implicit_tol=1e-10 (p1110/p12c never ran a link) | 6988.1 -> 7209.9 (221.8) | 19,750 | 22.71 | 4.14 | 970 | 8 | 5285 | 2251 | 0 | 22.8 % | rc=0 wall limit. **Corner jet onset 7115 d, jet at 1.5 R0 7160 d, 1.8 R0 7200 d** (RSG_BLOB_1113.md) |
| V2a | viper 12796543 | athena_he_gpu72_p_bc35bfca_p13c (= p12r + opt-in floor cap, bitwise p12r with keys absent) | binary only | 7210.5 -> 7390.5 (180.0) | 16,450 | 14.81 | 3.24 | 945 | 0 | 2336 | 734 | 0 | n/a | rc=137 STOP (tcpause) |
| V2b | viper 12796543 | p13c | none | 7390.5 -> 7395.7 (5.2) | 500 | 0.47 | 3.36 | 897 | 0 | 68 | 30 | 0 | n/a | rc=137 STOP (tcpause) |
| V2c | viper 12796543 | p13c | none | 7395.2 -> 7401.2 (6.0) | 600 | 0.59 | 3.54 | 867 | 1 | 87 | 37 | 0 | n/a | rc=137 STOP (tcpause) |
| V2d | viper 12796543 | p13c | none | 7400.2 -> 7425.3 (25.1) | 2,500 | 2.32 | 3.34 | 868 | 0 | 361 | 133 | 0 | n/a | rc=137 STOP (user pause) |
| R3 | Raven 32210941 (pHG holder) | athena_he_a100_nofma_t_d3b3f2e5_p1115 (prod-1115) | + time/dt_retry=true, **time2_stage_local_be=true** | 7425.3 -> 7469.2 (43.9) | 4,050 | 2.15 | 1.91 | 937 | 2 | 739 | 246 | 0 | n/a | our scancel 11-13 01:16 PST (user stop: corner jet) |
| R4a | Raven 32089203 | athena_he_a100_nofma_a_1d197280_hgdtd6 (prod-1118) | rerun from rsg.01451 (7100.2 d, pre-jet-at-1.5R0); + D6/face_fct/face_mmat/zalesak cap4/lat_posfix=lag/lzx/eos_h2/tfloor 500 ..., k1_redist 0, local BE off | 7100.2 -> 7100.2 | 20 | 0.1 | - | dt 111 | 8 | 7 | 5 | 0 | n/a | killed after 6 min (05_rsgfix_FIX), relaunched with identical keys |
| R4b | Raven 32089203 | prod-1118 1d197280 | identical to R4a | 7100.2 -> 7156.4 (56.2) | 6,700 | 6.08 | 3.26 | 725 | **334** | 1082 | 527 | 0 | n/a | rc=137 killed (cost regression, RSG_SLOW/RSG_CAUSE_1113.md) |

Not production (excluded): smokes, the 11-13 real_d6/poskeys/rsgcause/rsgslow test arms, and prod-1116 (only the 20-cycle energy audit on rsg.01517).

## 2. Per binary (rates per 1000 cycles)

| binary | machine(s) | sim days | cycles | NC /1k | NOT POS /1k | NA (BE redo) /1k | stops caused by the code | d sim per wall h | known physics problem |
|---|---|---|---|---|---|---|---|---|---|
| Caltech 23a55658 (guard 0) | Caltech | 1821.6 | 167,443 | n/a (logs not here) | n/a | n/a | none recorded in runH.cfg/inalloc.out | n/a | none recorded (no corner/vertex note after 1675 d) |
| Caltech earlier prodH binaries (e04fbd50 .. 6c76e282) | Caltech | 2830 | 258,855 | n/a | n/a | n/a | none recorded | n/a | none recorded after the vtxfix restart |
| 698b8241 prod-1107v | Raven | 123.9 | 10,750 | 0 | 157 | 44.5 | none (rc 143 = external) | 26.2 | none seen |
| 864d120d prod-1109 | Raven | 542.0 | 48,000 | 0.02 | 196 | 65.1 | none | 24.1 | none seen (jet onset is 127 d after its last cycle) |
| 4e9f602c p12r | viper | 221.8 | 19,750 | 0.41 | 268 | 114.0 | none | 9.8 (MI300A) | **cube-corner jet launched in this link** (onset 7115 d) |
| bc35bfca p13c (bitwise p12r) | viper | 214.8 | 20,050 | 0.05 | 142 | 46.6 | none (all ends = pauses) | 11.8 (MI300A) | ran on the post-jet state; jet kept growing to 2.44 R0 |
| p12r + p13c together | viper | 436.6 | 39,800 | 0.23 | 205 | 80.0 | none | - | corner jet |
| d3b3f2e5 prod-1115 (+ local BE) | Raven | 43.9 | 4,050 | 0.49 | 182 | 60.7 | none (our scancel) | 20.4 | prod-1115 energy residual 3.6e-5 OK; local BE on prod-1116 CREATES energy (P1116_ENERGY_1113.md) |
| 1d197280 prod-1118 | Raven | 56.2 | 6,720 | **50.9** | 162 | 79.2 | killed for cost (2.1x s/cyc of prod-1109 on the same machine, mean dt 725 s) | 9.2 | Picard limit cycle across e_floor in a lateral pocket (RSG_CAUSE_1113.md: lat_posfix=lag, k1_redist 0, local BE off, face_mmat) |

Zero FATAL in every production log. No segment ended because of the code except prod-1118 (cancelled for its cost); all other ends are wall limits, pauses or our own cancels.

## 3. Conclusions

1. Most stable measured RSG binary: **prod-1109 864d120d** (Raven run.cfg.32026500): 542 d / 48k cycles in one link, NC 1, the lowest BE share measured (13 %), clean wall-limit exit, 24 d per wall hour on 6 Raven nodes. prod-1107v (698b8241, same key set minus 2 hydro keys) is equally clean but only 124 d.
2. Second: **p13c bc35bfca** (= p12r bitwise with the floor-cap keys absent; viper run.cfg.12796543): 215 d, NC 0.05/1k, NA 47/1k, all ends were pauses. Its p12r twin produced the corner jet, so p12r/p13c with the geo2/EMA128 keys is not cleared.
3. Not for pRSG: **prod-1118 1d197280** (NC 51/1k, 2x slower, the cause is in its new D6 keys), **prod-1116** (local BE creates energy), and **any binary with time2_stage_local_be=true** until the energy fix. prod-1115 has only 44 d of RSG evidence (with local BE on), too short to rank.
4. Recommended pRSG = the R2 key set: binary athena_he_a100_nofma_o_864d120d_fmadoff_p1109 (Raven) / athena_he_gpu72_q_864d120d_p1109 (viper), the KEYS lines of /raven/ptmp/jinma/raven_gap_1105/rsgE/run.cfg.32026500 (the r8 line + `vet_gd_tangential=true` + `mesh/layout_invariant=true` + `hydro/wall_wb_onesided=true` + `hydro/split_c2p_once=true`), with **no** geo2/EMA128/rotate keys, no `time2_stage_local_be`, no `time/dt_retry=true`. If a newer binary is wanted (rstabe restart state, speed), use p13c bc35bfca or prod-1115 d3b3f2e5 with exactly this key set and `time2_stage_local_be=false`, and gate it first by a 20-cycle comparison against 864d120d on rsg.01451 (not run here).

Caveats:
- "Stable" for prod-1109 may be the era, not the binary: 6445-6988 d had no ejection, while p12r ran into a real ejection (7130-7300 d, 3e29 g) plus the jet; the NOT POS / NA rates doubled there (268/114 vs 196/65), which can come from the state or from the keys.
- The jet link (V1) is the first link with geo2 + EMA128 + 3-D rotation and implicit_tol 1e-10 (key change at 6988 d, onset 7115 d). RSG_BLOB_1113.md says the jet is "not tied to any machine, binary or key change", but it only checked the binary switch at 7210 d and the seam keys; the rotation-key change 127 d before onset is untested. Hypothesis only.
- Seam guard: Caltech ran 1822 d with `cs_seam_rho_guard=0` and no recorded corner event after 1675 d; all our links since 11-07 use guard 4 (jet 790 d after the switch). Untested as a cause; the Caltech NC/NP counts are not available here (logs stayed on Caltech).
- RSG_CAUSE_1113.md calls run.1111_085731.log a "prod-1115" log; cmd.12796542.txt shows it was p12r 4e9f602c (prod-1115 had not been built yet).
- s/cycle across machines is not comparable (A100 x24 vs MI300A x8); within Raven: prod-1107v 1.58, prod-1109 1.69, prod-1115 1.91 (with local BE), prod-1118 3.26.
