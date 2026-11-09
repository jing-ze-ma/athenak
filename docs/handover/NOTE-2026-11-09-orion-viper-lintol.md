# NOTE 2026-11-09 (Orion -> viper): lintol-1109, which Picard stop test removes the R136 accuracy loss

Reply to `TASK-2026-11-09-orion-lintol.md`.

## Short answer
1. **Reproduced** viper's table on this CPU build: the passes are identical and the errors agree to the 2-3 digits printed (sec. 2).
2. **Almost all of the small-deck D error comes from the fresh-start IC relaxation, not from the 8 cycles.** The first
   41 implicit solves of a fresh run are the `ic_rad = fs_relax` pseudo-time steps (relax_n = 40, before cycle 0). They take
   104 of D's 152 passes and 26 of CE0's 28 extra passes. The relaxed IC (nlim = 0) is already 3.8e-7 / 7.5e-9 / 8.3e-7 off R2's;
   R3's relaxed IC is 1e-11 off. Restarted from one common IC (R2's relaxed state), the 8 cycles give **D = CE0 = every guard
   variant to every printed digit, at 1.8e-9 / 3.6e-9 / 2.5e-9**, which is R3 level. In a 108-cycle restart (t = 0 .. 1500 s)
   the final states of D, CE0, CEg1.5/2/3, CEnf, TRM3e-5/1e-5/3e-6 are also identical to every printed digit.
3. **Why the stop test does not matter in these cycles:** `implicit_cs2_vimp_npass = 3` (conservative2 + vimp, i.e. every
   cycle solve) forces 3 passes. The 3rd-pass change is ~1e-5, 10x below tol_eff = 1e-4, so the stop is decided by the plain
   test, never by the contraction estimate. Only the 2-3 cold-start solves after the restart reach a 4th/5th pass, which is
   where CE0's +2 passes come from. **This deck therefore cannot rank stop tests for rst8 cycles.** On rst8 the estimate must
   be stopping solves at pass >= 4, which never happens here.
4. In the phase where the estimate does stop (the relaxation solves), **no existing-key setting cheaper than CE0 meets the
   1.5e-8 target**: guard 1.5 / 2 / 3 give 9e-8 / 1.9e-7 / 5e-7 in the relaxed IC, saving only 3 / 6 / 16 of CE0's 26 extra passes; nofloor changes nothing; tol_rel_max 1e-5 gives 8e-8; tolcap and an exact
   inner solve (ew_max 0) do not help. The cheapest setting that meets the target is `implicit_conv_est_guard = 1` (or `<= 1.2`),
   which is CE0 here.
5. **Long runs: the residual error is the inexact inner solve, not the Picard stop** (sec. 5).

## 1. Build, runs
- Binary: the existing MPI CPU Release he_star_m1 build of exactly **f34ada236** (`git rev-parse HEAD` checked) from the
  r136fast task, `/orion/ptmp/jinma/r136fast_1109/wt_base/build_base/src/athena`. It was built with
  `CMAKE_CXX_FLAGS=-ffp-contract=off` (gcc/13, openmpi/4.1, Kokkos SERIAL, Release -O3). It was copied, not rebuilt.
- `run_lintol.sh` as is, except the launcher: `srun -n 4 --exact --mem=12G --cpu-bind=cores` inside a p.shared
  `--qos=serial` allocation, with all arms in parallel on one node (mpirun is not on Orion's PATH). The deck copies come from
  the same `mk()`, with keys REPLACED. Each arm takes ~2 min (fresh) or ~1 min (8-cycle restart) on 4 cores, and ~6-7 min for 108 cycles.
- Work dir: `/orion/ptmp/jinma/lintol_1109` (driver `scripts/drive.sh`, arm lists `scripts/arms_*.txt`, ACC.txt per set;
  bins and rst deleted after acc.py; base/ACC.txt was re-typed from the acc.py console output after a slip deleted the file).
- Sets: **base/x1/x2** = fresh, tlim 120 (as the TASK); **xic** = fresh, nlim 0 (relaxed IC only); **xr** = restart from
  R2's relaxed IC (`ic/R2ic/rst/r136.00001.rst`, t = 0), tlim 120 (8 cycles, 15 solves); **xL** = the same restart, tlim
  1500 (108 cycles, 226 solves). The restart warns "no predictor increment / no hesdirk2 slope": the first step cold starts,
  identically in every arm.
- Error = acc.py E max per region vs R2 (fresh, xic, xr). For xL the reference is R3, because R2 vs R3 is itself 4e-8 there.
- Passes / inner = `bicgstab: outer passes=... inner iterations mean=... total=...`. "vs D" = extra outer passes.

## 2. Reproduction (fresh, tlim 120, vs R2)
| arm | interior | transition | thin | passes | viper |
|---|---|---|---|---|---|
| R3 | 3.0e-9 | 1.5e-9 | 1.4e-9 | 224 | 3.0e-9 / 1.5e-9 / 1.4e-9 |
| D | 4.1e-7 | 3.3e-8 | 2.8e-7 | 152 | 4.1e-7 / 3.3e-8 / 2.8e-7, 152 |
| Kall | 4.1e-7 | 3.2e-8 | 2.8e-7 | 151 | 4.1e-7 / 3.2e-8 / 2.8e-7, 151 |
| CE0 | 1.4e-8 | 3.6e-9 | 5.8e-9 | 180 | 1.4e-8 / 3.6e-9 / 5.8e-9, 180 |
| CEg1 (viper) / CEg1.2 (here) | = CE0 | = CE0 | = CE0 | 180 | = CE0, 180 |

All match.

## 3. Fresh, tlim 120: all arms (vs R2; E max per region; outer passes; inner mean / total)
| arm | interior | transition | thin | outer passes | vs D | inner mean / total |
|---|---|---|---|---|---|---|
| R3 | 3.0e-09 | 1.5e-09 | 1.4e-09 | 224 | +47.4 % | 3.64 / 816 |
| D | 4.1e-07 | 3.3e-08 | 2.8e-07 | 152 | +0.0 % | 1.37 / 208 |
| Kall | 4.1e-07 | 3.2e-08 | 2.8e-07 | 151 | -0.7 % | 1.64 / 248 |
| CE0 | 1.4e-08 | 3.6e-09 | 5.8e-09 | 180 | +18.4 % | 1.48 / 266 |
| CEg3 | 6.5e-07 | 2.1e-08 | 1.5e-07 | 164 | +7.9 % | 1.41 / 232 |
| CEg10 | 4.7e-07 | 3.2e-08 | 2.0e-07 | 155 | +2.0 % | 1.37 / 213 |
| CEnf | 4.1e-07 | 3.3e-08 | 2.8e-07 | 152 | +0.0 % | 1.37 / 208 |
| CEg1p2 | 1.4e-08 | 3.6e-09 | 5.8e-09 | 180 | +18.4 % | 1.48 / 266 |
| CEg1p5 | 1.1e-07 | 5.6e-09 | 1.9e-08 | 177 | +16.4 % | 1.46 / 259 |
| CEg2 | 2.4e-07 | 1.2e-08 | 4.9e-08 | 174 | +14.5 % | 1.45 / 253 |
| CEnfg3 | 6.5e-07 | 2.1e-08 | 1.5e-07 | 164 | +7.9 % | 1.41 / 232 |
| TR03 | 3.9e-07 | 3.1e-08 | 9.9e-08 | 158 | +3.9 % | 1.39 / 220 |
| TR01 | 3.9e-07 | 2.9e-08 | 9.5e-08 | 159 | +4.6 % | 1.40 / 222 |
| TRM3e5 | 5.1e-07 | 2.5e-08 | 1.2e-07 | 160 | +5.3 % | 1.41 / 226 |
| TRM1e5 | 9.7e-08 | 6.3e-09 | 2.1e-08 | 175 | +15.1 % | 1.47 / 258 |
| EW1e3 | 4.1e-07 | 3.1e-08 | 2.8e-07 | 152 | +0.0 % | 1.59 / 242 |
| EW0 | 4.1e-07 | 3.1e-08 | 2.8e-07 | 151 | -0.7 % | 4.45 / 672 |
| FRSr1 | 2.9e-07 | 1.0e-07 | 4.9e-08 | 174 | +14.5 % | 1.45 / 252 |
| CE0TR2 | 1.4e-08 | 3.6e-09 | 5.8e-09 | 180 | +18.4 % | 1.48 / 266 |
| TCL1 | 4.1e-07 | 3.3e-08 | 2.8e-07 | 152 | +0.0 % | 1.37 / 208 |
| R2 | (ref) | (ref) | (ref) | 195 | +28.3 % | 4.04 / 787 |

Notes. TR03/TR01 (implicit_tol_rel 0.03 / 0.01) and CE0TR2 (tol_rel 0.2) barely move, because tol_eff = min(tol_rel r0,
tol_rel_max) is capped at 1e-4 in almost every solve. TCL1 = implicit_lin_tolcap_last 1 is bitwise D. EW0 (exact inner solve,
D stop) is D's error at 3x the inner iterations, so the loss is not inner inexactness. FRSr1 (implicit_final_resid_stop_rel 1,
not recommended in picstop-1017) makes the error worse.

Opacity-lag check (x2; each D/CE0 vs an R2 with the same lag):
| lag | D | CE0 | passes D / CE0 / R2 |
|---|---|---|---|
| implicit_opac_lag_passes 1 | 2.1e-12 / 7.3e-12 / 8.3e-12 | - | 132 / - / 129 |
| 2 | 4.0e-8 / 1.5e-8 / 7.2e-8 | 4.0e-8 / 1.5e-8 / 7.2e-8 | 148 / 174 / 169 |
| 3 (deck) | 4.1e-7 / 3.3e-8 / 2.8e-7 | 1.4e-8 / 3.6e-9 / 5.8e-9 | 152 / 180 / 195 |
| none (opacities every pass) | 4.1e-7 / 3.6e-8 / 2.6e-7 | 4.5e-8 / 7.5e-9 / 2.4e-8 | 161 / 192 / 210 |

So the loss in the relaxation solves follows the opacity nonlinearity. With opacities frozen after pass 0 (lag 1) any stop is
exact. With lag 3, a 2-pass stop (the estimate's typical stop in relax solves 6-27: changes 2.7e-2 -> 1.6e-3, q = 0.06,
estimate 1.0e-4 just below tol_eff) ends before the pass whose opacities are frozen, so it is a different discretisation from the
reference's. The estimate itself is accurate there: CE0's 3rd-pass change is 8.5e-5, as predicted.

## 4. Where the stop test binds: the relaxed IC only (xic, fresh nlim 0, vs R2)
| arm | interior | transition | thin | outer passes | vs D | inner mean / total |
|---|---|---|---|---|---|---|
| R3 | 1.9e-12 | 7.8e-12 | 1.7e-11 | 164 | +57.7 % | 3.78 / 620 |
| D | 3.8e-07 | 7.5e-09 | 8.3e-07 | 104 | +0.0 % | 1.35 / 140 |
| CE0 | 2.0e-08 | 6.4e-09 | 7.2e-09 | 130 | +25.0 % | 1.52 / 197 |
| CEg1p5 | 9.2e-08 | 2.5e-08 | 1.5e-07 | 127 | +22.1 % | 1.50 / 190 |
| CEg2 | 1.9e-07 | 5.5e-08 | 3.5e-07 | 124 | +19.2 % | 1.48 / 184 |
| CEg3 | 5.0e-07 | 1.0e-07 | 5.7e-07 | 114 | +9.6 % | 1.43 / 163 |
| TRM1e5 | 8.0e-08 | 2.5e-08 | 1.5e-07 | 127 | +22.1 % | 1.50 / 190 |
| EW0 | 3.8e-07 | 3.8e-09 | 8.6e-07 | 103 | -1.0 % | 4.74 / 488 |
| R2 | (ref) | (ref) | (ref) | 147 | +41.3 % | 4.10 / 603 |

Here CE0 is the only setting at the target. Guard 1.5 / 2 / 3 save 3 / 6 / 16 of CE0's 26 extra passes and are 5x / 10x / 25x
above the target. Tightening tol_rel_max to 1e-5 with the estimate on costs as much as guard 1.5 and is no more accurate.

## 5. Cycles only: restart from one common IC
8 cycles (xr, vs R2):
| arm | interior | transition | thin | outer passes | vs D | inner mean / total |
|---|---|---|---|---|---|---|
| R3 | 3.0e-09 | 1.5e-09 | 1.4e-09 | 60 | +25.0 % | 3.27 / 196 |
| D | 1.8e-09 | 3.6e-09 | 2.5e-09 | 48 | +0.0 % | 1.42 / 68 |
| Kall | 1.4e-09 | 3.1e-09 | 1.1e-09 | 48 | +0.0 % | 1.81 / 87 |
| CE0 | 1.8e-09 | 3.6e-09 | 2.5e-09 | 50 | +4.2 % | 1.38 / 69 |
| CEg1p5 | 1.8e-09 | 3.6e-09 | 2.5e-09 | 50 | +4.2 % | 1.38 / 69 |
| CEg2 | 1.8e-09 | 3.6e-09 | 2.5e-09 | 50 | +4.2 % | 1.38 / 69 |
| CEg3 | 1.8e-09 | 3.6e-09 | 2.5e-09 | 50 | +4.2 % | 1.38 / 69 |
| CEg10 | 1.8e-09 | 3.6e-09 | 2.5e-09 | 49 | +2.1 % | 1.41 / 69 |
| CEnf | 1.8e-09 | 3.6e-09 | 2.5e-09 | 48 | +0.0 % | 1.42 / 68 |
| TRM3e5 | 1.8e-09 | 3.6e-09 | 2.5e-09 | 48 | +0.0 % | 1.42 / 68 |
| TRM1e5 | 1.8e-09 | 3.6e-09 | 2.5e-09 | 48 | +0.0 % | 1.42 / 68 |
| TR01 | 1.8e-09 | 3.6e-09 | 2.5e-09 | 48 | +0.0 % | 1.42 / 68 |
| EW1e3 | 1.3e-09 | 3.0e-09 | 2.2e-09 | 48 | +0.0 % | 1.56 / 75 |
| FRSr1 | 2.4e-07 | 9.6e-08 | 4.4e-08 | 48 | +0.0 % | 1.42 / 68 |
| R2 | (ref) | (ref) | (ref) | 48 | +0.0 % | 3.83 / 184 |

108 cycles (xL, vs **R3**):
| arm | interior | transition | thin | outer passes | vs D | inner mean / total |
|---|---|---|---|---|---|---|
| R3 | (ref) | (ref) | (ref) | 903 | +32.0 % | 3.22 / 2912 |
| R2 (vs R3) | 3.9e-08 | 1.5e-08 | 5.1e-08 | 793 | +15.9 % | 3.49 / 2770 |
| D | 6.9e-08 | 5.1e-08 | 5.5e-08 | 684 | +0.0 % | 1.34 / 914 |
| CE0 | 6.9e-08 | 5.1e-08 | 5.5e-08 | 686 | +0.3 % | 1.33 / 915 |
| Kall | 6.9e-08 | 5.0e-08 | 5.5e-08 | 684 | +0.0 % | 1.38 / 943 |
| EW1e3 | 6.9e-08 | 5.0e-08 | 4.4e-08 | 684 | +0.0 % | 1.37 / 934 |
| EW1e4 | 6.9e-08 | 2.9e-08 | 5.3e-08 | 684 | +0.0 % | 1.97 / 1349 |
| EW0 | 6.9e-08 | 2.9e-08 | 5.4e-08 | 684 | +0.0 % | 3.75 / 2563 |
| TRM1e5 | 6.9e-08 | 5.1e-08 | 5.5e-08 | 684 | +0.0 % | 1.34 / 914 |
| CE0_TRM1e5 | 3.4e-09 | 4.7e-08 | 4.6e-08 | 902 | +31.9 % | 1.40 / 1267 |
| TRM1e6 | 3.9e-08 | 4.8e-08 | 5.3e-08 | 795 | +16.2 % | 1.42 / 1127 |
| CE0_TRM1e6 | 3.0e-09 | 4.6e-08 | 5.2e-08 | 906 | +32.5 % | 1.40 / 1271 |
| EW0_TRM1e6 | 3.9e-08 | 1.5e-08 | 5.1e-08 | 793 | +15.9 % | 3.49 / 2770 |
| EW0_CE0_TRM1e5 | 5.1e-10 | 4.6e-10 | 2.7e-08 | 901 | +31.7 % | 3.23 / 2910 |

- In the 8- and 108-cycle restarts every stop-test variant is identical to D to every printed digit (CE0 +2 passes =
  +0.3 % over 108 cycles; the 2 passes are in the cold-start solves).
- Over 108 cycles the interior error (7e-8) goes only with a Picard tolerance below the 3rd-pass change, i.e. a 4th pass in
  nearly every solve: CE0 + tol_rel_max 1e-5 gives 3e-9 at +32 % passes.
- The transition/thin error (~5e-8) is unchanged by any Picard setting, even at R3's pass count (CE0_TRM1e6 906 passes vs R3 903).
  It goes with the inner tolerance: EW1e4 (implicit_lin_ew_max 1e-4) and EW0 bring transition from 5.0e-8 to 2.9e-8, and EW0 + Picard tol_rel_max 1e-6 (= R2 bitwise) to 1.5e-8. EW0 + CE0 + tol_rel_max 1e-5 gives 5e-10 / 5e-10 / 2.7e-8. The thin maximum (~5e-8) stays at R2-vs-R3 level, i.e. at the uncertainty of the R2 reference itself after 108 cycles. So in long runs the 1.5e-8 target is a Picard + inner TOLERANCE question, not a stop-test question, and is below what R2 resolves; the stop variants do not touch it.

## 6. Code-change ideas (described, not committed, not tested)
- The relax-phase failure is a 2-pass stop whose q comes from passes (0,1). Pass 0's change is the predictor error, not a
  contraction of the fixed-point map, and under implicit_opac_lag = pred with N = 3 passes the map changes after pass 2. A minimal
  variant: allow the estimate stop only for `it + 1 >= impl_opac_lagn` when `impl_opac_lagp` (one line next to the guard at
  rad_m1_implicit.cpp ~16654, `if (pc && impl_opac_lagp && it + 1 < impl_opac_lagn) pc = false;`). In the cycles here
  (npass 3) it costs nothing, and in the relax solves it equals CE0 in solves 6-27. It does nothing for the pass >= 4 stops that rst8 must have.
- The contraction requirement q < 0.5 (code) could be tightened, or the estimate given a safety factor
  (`resid*q/(1-q)*S < tol_eff`), like `time2_one_pass_safety` for the one-pass stop. There is no existing key for either.
  In the relax solves q ~ 0.05, so a tighter q bound would not catch them, and S ~ 1.2 would.

## 7. Recommendation
- On this deck, the only setting that meets <= 1.5e-8 wherever the stop test binds is **`implicit_conv_est_guard = 1`**
  (bitwise `implicit_conv_est = 0`; guard 1.2 is also bitwise CE0 here):
  ```
  <rad_m1>
  implicit_conv_est_guard = 1
  ```
  Its pass overhead over D: +18.4 % fresh (152 -> 180), but **+26 of those +28 passes are IC-relaxation passes**. In
  evolution cycles: +4.2 % over 8 cycles (48 -> 50) and **+0.3 % over 108 cycles (684 -> 686)**. The cycle overhead here is
  all from the cold-start solves.
- No cheaper existing-key setting exists on this evidence. Guard >= 1.5, nofloor, tolcap (all/last), tol_rel, tol_rel_max
  >= 1e-5 and ew_max all fail the target where the estimate binds, and are identical to D where it does not.
- Viper's rst8 +5.2 % for CE0 must come from solves that need > 3 passes, which this small fresh deck never reaches in its
  cycles. To rank cheaper variants on CPU, an rst8-like state is needed (e.g. a 1-D column or a small sp sector of the rst8
  state). With one, guard 1.5 / 2 and the opac-lag line above are the candidates to test first.
