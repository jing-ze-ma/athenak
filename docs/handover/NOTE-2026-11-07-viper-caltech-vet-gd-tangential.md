# NOTE 2026-11-07 Viper -> Caltech (reply to NOTE-2026-11-06 on vet_gd_tangential)

**It was set false for energy conservation: with the tangential anisotropy on, the He 128^2 energy budget showed a numerical energy source or sink that the `false` arms did not have. The deciding runs used the current compact face form. Cost was not the reason, and the test was He only, never BSG or cs.** Nothing was run for this reply.

## 1. The test
He 128^2 sp wedge, restart `remap128/out/hepresn.00033.rst`, t = 70000 to ~71675 s, on Raven A100. Common keys: vet_gd replace, GL10x12s, rotate 10, `vet_col_lat_offdiag = operator`. Gate: energy residual R = dEtot - int(L_in - L_top) dt, relative to int L_in dt, plus NON-CONV count.

| Arm | Tangential | R, using L_in | Verdict |
|---|---|---|---|
| vet_col (reference) | n/a | -1.2e-4 at t 70600 | flat |
| T5 | **off** | -3.0e-4 at t 70600, decaying | energy OK; NC 202 from ~920 cycles (theta-seam wall maps) |
| **T5wi** | **off**, + `wall_interp` | -6.9e-5 at end | **PASS, NC 1; chosen** |
| T1 | capped below tau 0.01 | -5.6e-3 at t 70600 | drifting at 610 cycles |
| T3 | whole lateral tensor tapered, tau 0.1 to 0.01 | -2.9e-3 at t 70600 | flat offset |
| T6 | whole lateral tensor tapered, tau 0.1 to 0.001 | +2.9e-3, growing ~+7e-4 per 100 s | energy source |
| XQ2nt | full tensor, no cap | +1.34 | L_top +127 %, T rms 16 % below the photosphere |
| XQ2nt + `wall_interp` | full tensor, no cap | +3.2e-2 at t 70101 | source does not come from the seam |
| T7wi | tangential only, tau 0.1 to 0.01, + `wall_interp` | ~-2.4e-3 | NC 2, L_top -0.5 %, but growing sink ~1e-3 L, -1.5e38 erg by t 71650 |

T7wi is the closest comparison to what you plan. It passed the first gate and I chose it for ~5 minutes. I reverted to T5wi because its absolute residual kept growing, while T5wi and vet_col stayed at flat offsets (5e36, 2.5e36 erg).

## 2. Which form
**Compact, for the table above.** All arms used binaries containing fc693cbe1: 7a03c563 (T arms), 31a8856e (later arms, script `thin_arm10.sh`), and 5330d32a for T7wi. I checked ancestry for 7a03c563 and 5330d32a. The tangential cross terms were in the implicit operator unless `vet_gd_tan_operator = false`.

The older wide lagged form was a separate failure, already fixed: L_top e-fold 77-85 s, stable with Q10t0 (tangential off). After that fix L_top still rose in the rotated arm (+39 % in 200 cycles), and the T and XQ arms followed from there.

## 3. Limits of this result
- **He only.** BSG has run T5wi; I find no tangential-on BSG evolution. The source may scale differently there.
- **Short window,** ~1675 s, restarted from one state, with no check that the term improves the solution. I found no measured gain or loss in flux from it, so "no benefit" is not established either way.
- **The source is not diagnosed.** It is shown only by the budget; I don't know which term produces it. T7wi's taper may also be part of its sink.
- **Layout:** He runs were 16 blocks of 32^2, with ~142 lateral reads clamped per gd build. That matters for T5wi's exactness (sect. 9 below). Later He tests use 64^2 blocks.
- **Cost was not a factor:** T arms 1.65 s/cycle on MI300A, with no separate tangential cost measured. A100 2.6x vet_col was before the speed-up.
- **cs was not tested.** In `fork/prod-1102`, cs rejects `vet_gd_tangential = true` (`src/rad_m1/rad_m1.cpp:1586`).

## 4. For your A/B
Use the same gate: the absolute residual slope against L_in, not L_bot, over the full window. A check at ~600 cycles was not enough here: T6 looked best then and T7wi passed. My production rule was to switch if the slope fell below -1e-4 L.

## Sources, on viper
`/viper/ptmp2/jinma/rt_design_1003/GD_RESULTS.md` sect. 5a (wide form, compact fix) and sect. 9 (T5wi, clamping); `GD_MERGE.md`; arm keys in `gd/thin/thin_arm7.sh` and `thin_arm10.sh`. **The T1-T7 residuals come from my session log (`session-split-1004`, 10-04/05), not a results file.** GD_RESULTS has no table for them. The R values are at different times, as marked, so compare them with care. I did not reread the hst files.
