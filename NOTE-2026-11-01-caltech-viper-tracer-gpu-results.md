# NOTE 2026-11-01 Caltech -> viper (USER): tracer GPU tests T0-T2 PASS; T1 xe arm cuts RKG stages 17 -> 7 at 0.4 rot; floor tracer nets ~36 % of the floor mass

Reply to TASK bf54375e6 and NOTEs dc10d2cc / 015b0abd / 6877cc88.

## Binaries and runs
All builds are H200, fmad off (`-ffp-contract=off`, nvcc `--fmad false`), CUDA 12.9, hpcx 2.17.1, incremental build in builds_rstfixx/inc_gpu_dhj_fmadoff_mhdfloor (the tree that built b935e145).

| tag | commit | binary (builds_rstfixx/bin/) | md5 |
|---|---|---|---|
| ref | fin51-r-fix b935e145 | athena_gpu_dhj_fin51_b935e145_fmadoff | 1b7f60a0d65a6fa16553f5c86341ba41 |
| test | tracers-test-1101 2b72d789 | athena_gpu_dhj_trc_2b72d789_fmadoff | bb954247b6671b9e218e58f1b4e17e90 |
| test2 | tracers-test-1101 89213fcb | athena_gpu_dhj_trc_89213fcb_fmadoff | 16b9fb5be8069c8e778852da94ec7bf4 |
| def | tracers-defaults-1101 ac98652f | athena_gpu_dhj_trcdef_ac98652f_fmadoff | 60cc2613855548a12651007675c868b8 |

- Every run uses 1 rank on 1 H200, with BVALS_CC_SPLIT=1. Run dirs: /resnick/groups/carnegie_poc/jingze/trc_1101/runs/.
- **Deck:** trc_1101/decks/r1m_trc.athinput. It is mhdfloor_1101/r1m_prod.athinput (top_bc outflow1, rkg_dt_safety 0.5, rkg_lean true, flux_hst_floor true) plus:
  - cs_seam_flux = shared;
  - am_hst_split = false;
  - valfven_noinfall_fac = 4;
  - explicit `xe_tracer = false` and `floor_tracer = false` in `<mhd>`, which the arms override.
- **Arms run concurrently,** 2-3 per GPU at first. Wall times are therefore only indicative (see T1 cost).

## NHISTORY check (your "PLEASE CHECK r1m TODAY")
- **YES, the r1m production deck has the combination,** and it crashes.
  - r1m_prod.athinput has am_hst_split = true, flux_hst_floor = true and roche = true.
  - On b935e145 that makes 9 + 3 + 12 + Ttide = 25 user columns.
  - The deck **segfaults at initialisation** in DhjFluxHistory (std::string assign, label[24]) from HistoryOutput::LoadOutputData in Driver::Initialize, before cycle 0 (run dir t0a_ref_split_segv).
- **No production run was hit:** r1m (w121r1_1101) was stopped at 00:22, on an older binary without Ttide.
- **Every relaunch on b935e145 or later needs am_hst_split = false.** The defaults branch drops Mdrain/Edrain with a NOTE, but the Ttide overflow is in the base code.

## T0. Key off bitwise: PASS

### (a) r1m deck from t = 0
- Deck: r1m_prod with am_hst_split = false, because of the segfault above. tlim 2000 s = 186 cycles; rst and bin at 0, 1000 and 2000 s.

| comparison | rst payload (3 files) | dhj.mhd.hst | bin data (3 files) | dhj.log | dhj.user.hst |
|---|---|---|---|---|---|
| b935e145 vs 2b72d789 (keys absent) | identical (whole files identical) | identical | identical | identical | Mfloor/Efloor/Efloor_rt differ <= 3.3e-15 rel; all other columns identical |
| b935e145 vs b935e145 (repeat) | - | - | - | - | same three columns differ <= 4.1e-15 rel (GPU atomics are base behaviour) |
| b935e145 vs ac98652f (explicit xe_tracer/floor_tracer = false) | identical payload | identical | identical data (the bin header carries the extra keys) | - | same three columns <= 5.5e-15 rel |

### (b) Restart check
- No developed W121 cs MHD restart was at hand. T1's off arm instead serves as the restart key-off check, and it is stronger:
  - **2b72d789/89213fcb off arm vs p_f4,** which ran on 81d407dc with the same keys, from t = 0 to 0.5 rot, 8774 cycles.
  - **Result:** dhj.mhd.hst identical, and all 5 bin dumps identical in data.
  - The counters, the floor columns and E_B agree to every printed digit. user.hst differs only because p_f4 had am_hst_split on.
- The csvertex A arm (9f9a3d50, key off) restarted from that file at 0.25 rot. Its dt reproduces p_f4/off to every printed digit (5.033 s at 0.4 rot, 4.507 s at 0.5 rot).

## T1. r1m production arms, t = 0 to 0.5 rot: PASS for robustness and ledger; floor tracer partial; xe saves stages

- **Binaries:**
  - off and xe: 2b72d789. Mdrain/Edrain do not exist in these arms, and the only change in 89213fcb is those columns.
  - fl and both: 89213fcb. 2b72d789 fl/both runs were stopped and restarted when the column fix arrived.
- **Data:**
  - hst values are means over 0.4-0.5 rot.
  - dt and RKG stages come from the run.log ndiag lines at the nearest 100-cycle mark.
  - Bands at 0.5 rot come from mhdfloor_1101/bands.py.

| arm | dt 0.25 / 0.40 / 0.50 rot [s] | RKG stages at 0.25 / 0.40 / 0.50 rot | dt_diff at 0.40 rot [s] | Mfloor net [g/s] | Mdrain [g/s] | gross = net + Mdrain | Mdot_top [g/s] | E_B 0.5 rot [erg] |
|---|---|---|---|---|---|---|---|---|
| p_f4 (81d407dc) | 6.028 / 5.033 / 4.524 | 17 / 16.7 / 16 | 0.0729 | 1.561e14 | - | 1.561e14 | 3.42e11 | 1.5284e32 |
| off | 6.028 / 5.033 / 4.524 (= p_f4 bitwise) | 17 / 16.7 / 16 | 0.0729 | 1.561e14 | - | 1.561e14 | 3.42e11 | 1.5284e32 |
| xe | 5.984 / 4.803 / 4.319 | 6.6 / 7.0 / 7.0 (range 5-11) | 0.495 | 1.432e14 | - | 1.432e14 | 3.44e11 | 1.6252e32 |
| fl | 6.026 / 5.036 / 4.520 | 17 / 16.7 / 16 | 0.0729 | 9.83e13 | 5.53e13 | 1.535e14 | 4.6e10 | 1.5276e32 |
| both | 5.983 / 4.800 / 4.318 | 6.5 / 6.7 / 7.2 (range 5-11) | 0.522 | 8.76e13 | 5.29e13 | 1.405e14 | 4.3e10 | 1.6235e32 |

### Robustness and counters
- **Robustness:** no FATAL, nan or NON-CONVERGED in any arm.
- **Counters** (sum over the run; eos_vceil, eos_fail, c2p_it and fofc are 0 in every arm):

| arm | eos_dfloor | eos_efloor | eos_tfloor |
|---|---|---|---|
| off | 2.553e9 | 2.506e8 | 3.242e8 |
| xe | 2.418e9 | 2.504e8 | 3.514e8 |
| fl | 2.607e9 | 2.521e8 | 3.151e8 |
| both | 2.474e9 | 2.423e8 | 3.218e8 |

- **Counter verdict: PASS.** Every arm is within +8 % of off; xe tfloor is the largest.

### Floor tracer range (0 <= s <= rho)
- **PASS.** The floor scalar (s/rho in the bin) lies in [0, 1.0000] in fl and both, with no cell < 0 or > 1.
- The x_e scalar lies in [0, 1.7e-4].

### Mass ledger: PASS (closes as in off)
- Test: dM(mhd.hst) vs the integral of (Mfloor - Mdot_top - Mdot_bot) dt, from t = 0 to 0.5 rot.

| arm | residual [g] | fraction of the integral of Mfloor | integral of Mdrain [g] |
|---|---|---|---|
| off | 5.4e14 | 7.3e-5 | - |
| xe | 5.5e14 | 7.6e-5 | - |
| fl | 5.2e14 | 9.5e-5 | 2.03e18 |
| both | 5.3e14 | 1.0e-4 | 2.02e18 |

- The residual is the same in every arm (hst sampling), so the drain booking closes.

### Net vs gross floor mass: PARTIAL
- Mdrain <= gross holds.
- **Net is not << gross:** net/gross = 0.64 (fl) and 0.62 (both). The tracer removes ~36-38 % of the floor-added mass at 0.4-0.5 rot. The per-row maximum Mdrain/gross over the run is 0.41 (fl) and 0.49 (both).
- **Side effect at the top:**
  - Mdot_top falls 8x (3.4e11 -> 4.5e10 g/s).
  - The mass-weighted v_r at p < 1e-6 bar goes from -148 m/s to -13 m/s.
  - M(p < 1e-6 bar) falls from 1.838e19 to 1.70e19 g.

### p > 1e-5 bar vs off: PASS from 1e-4 bar down; MARGINAL in the top shell
- **Measure:** horizontal mass-weighted shell means at 0.5 rot, maximum |rel diff| over shells, with p_f2 vs p_f4 as the fac 2/4 spread. The measured spread is larger than 0.1 % at the top.

| vs off | T, p > 1e-5 | rho, p > 1e-5 | T, p > 1e-4 | rho, p > 1e-4 | T, p > 1e-3 | rho, p > 1e-3 |
|---|---|---|---|---|---|---|
| spread p_f2 vs p_f4 | 2.1e-3 | 3.3e-3 | 3.7e-4 | 1.25e-3 | 1.7e-4 | 8.1e-4 |
| xe | 1.35e-3 | 1.05e-3 | 3.2e-4 | 4.6e-4 | 3.2e-4 | 4.6e-4 |
| fl | 4.9e-3 (at 1.2e-5 bar) | 6.4e-3 | 8.9e-4 | 1.3e-3 | 1.5e-4 | 2.3e-4 |
| both | 4.0e-3 (at 1.2e-5 bar) | 6.7e-3 | 8.8e-4 | 1.35e-3 | 2.2e-4 | 2.5e-4 |

- **xe** is inside the spread everywhere.
- **fl and both** are 2x the spread in the shells at 1.1-1.2e-5 bar, where the drain acts. They are at or near the spread at p > 1e-4 bar and inside it at p > 1e-3 bar.
- **B shell means** differ by 7 % (xe) and 12 % (fl) vs a spread of 9.5 %.

### Bands at 0.5 rot, 1e-5..1e-4 bar
p_f4 reference: 2.648e20 g / 2732 K / 2.53 km/s / 8.2 G.

| arm | mass [g] | T [K] | v_rms [km/s] | B_rms [G] |
|---|---|---|---|---|
| off | 2.648e20 | 2732 | 2.53 | 8.2 |
| xe | 2.664e20 | 2719 | 2.55 | 7.9 |
| fl | 2.646e20 | 2733 | 2.53 | 8.1 |
| both | 2.666e20 | 2717 | 2.55 | 7.9 |

### x_e payoff: RKG stage count falls 17 -> 7
- **Result:** the RKG stage count falls from 16.7 to 7.0 at 0.4 rot (-58 %). dt_diff grows from 0.073 s to 0.50 s, 6.8x, which means eta_max per dx^2 drops 6.8x.
- **This contradicts your 6877cc88 expectation of -3 to -6 %.** It is consistent in xe and both and over 0.25-0.5 rot.
- **The hydro dt is 4.6 % lower** (4.80 vs 5.03 s at 0.4 rot). E_B is 6 % higher.
- **Wall time:** off and xe ran side by side on one GPU from the same start. xe finished 0.5 rot in 25.6 min vs 35.1 min (-27 %; per cycle 0.170 vs 0.236 s). This is indicative only, because the GPU was shared.
- **Not done:** where the binding eta_max cell sits in off vs xe (no rkg_lmax_diag in these runs).

### Known defect (6877cc88) seen
- In xe at 0.5 rot, the night-side x_e at p < 1e-8 bar has median 2e-9, and 35 % of cells are < 1e-12. At 1e-8..1e-6 bar, 5 % are < 1e-12, with minimum 0.
- This is the floor-mass dilution you describe; the dayside is not affected (min 1e-6).

## T2. Restart subsets

### On 89213fcb, from T1 rst at 0.25 rot (cycle 3422/3426), +300 cycles, cycle-matched
| test | result |
|---|---|
| both-file, xe_tracer=false | clean, rc 0 ("x_e tracer is not read") |
| fl-file, xe_tracer=true | clean ("starts from local Saha") |
| off-file, both on | clean (Saha / zero starts) |
| both -> both, 300 vs 150 + 150 | **NOT bitwise**: payload differs from byte 232, dt 5.899234 vs 5.899354; log "no ck_sph_face = 5 VEF state in this file ... (this restart is not bitwise)". This is the known ck restart-slab bug (015b0abd). |
| control: off -> off, 300 vs 150 + 150 (89213fcb) | payload IDENTICAL (375138468 bytes after `<par_end>`) |

### On ac98652f (defaults branch): PASS
- Fresh base rst files come from t = 0, 600 cycles, with explicit keys.

| test | result |
|---|---|
| both-file, xe_tracer=false, to cycle 900 | clean, rc 0 |
| fl-file, xe_tracer=true, to 900 | clean (Saha start) |
| off-file, both on, to 900 | clean (Saha / zero starts) |
| off-file, floor_tracer on, nlim = 600 (0 cycles; the segfault case) | rc 0, rst written |
| both -> both, 300 vs 150 + 150 | **PAYLOAD IDENTICAL** (392956068 bytes after `<par_end>`); time, dt and ncycle identical (9395.80328382235, 10.1718365622474, 900) |

- The whole-file cmp differs only in the header (input text). I compared the payload after `<par_end>`, because the header length changes with the nlim text, so a fixed tail-20000 cmp misaligns.

## T3 (optional; partial at the time of writing)
| arm (0.5 rot, 89213fcb) | Mfloor net | Mdrain | gross | net/gross | dt 0.4 rot |
|---|---|---|---|---|---|
| fl default (margin 2, fmax 0.01) | 9.83e13 | 5.53e13 | 1.535e14 | 0.64 | 5.036 |
| fmax 0.03 | 9.84e13 | 5.54e13 | 1.538e14 | 0.64 | 5.041 |
| margin 4 | 1.163e14 | 3.59e13 | 1.521e14 | 0.76 | 5.032 |

- **fmax does not bind:** 0.01 and 0.03 give the same result.
- **margin 4 drains 35 % less.**
- The ledger closes in all of these arms (residual ~1e-4).
- fmax 0.003 and the xe alpha_fac 3 / 10 arms (1 rot) are still running. I will add them only if they change the picture.

## Summary
- T0: PASS (key off bitwise; floor columns only at GPU-atomic round-off, as in base).
- T1: PASS for robustness, counters, 0 <= s <= rho and the mass ledger.
  - Floor tracer: net/gross 0.62-0.64, so net is not << gross; Mdot_top falls 8x; the 1.2e-5 bar shells are 2x the fac 2/4 spread.
  - x_e tracer: 58 % fewer RKG stages, -4.6 % hydro dt, about -27 % wall time (indicative).
- T2: PASS on ac98652f. On 89213fcb the both-tracer continuation is not bitwise (the known bug); the no-tracer continuation is bitwise (single-tracer continuation not run on 89213fcb).
- NHISTORY: the r1m production deck segfaults on b935e145 at init; use am_hst_split = false.
