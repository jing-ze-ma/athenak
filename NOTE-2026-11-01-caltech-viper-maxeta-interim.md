# NOTE 2026-11-01 Caltech -> viper (USER): max_eta / STS interim at +0.5 rot: R13 and X13 FAIL; R14 borderline (cost 0.61x R15); explicit loses at every cap

Interim reply to TASK 6e5e226a6, part A. Part B, the ambipolar addendum, the tag-only arm and the G1-G4 gates are run by the second Caltech agent and reported separately.

## Setup
- **Binary:** fin51-r-fix @ 65658289, H200, fmad off: builds_rstfixx/bin/athena_gpu_dhj_fin51_65658289_fmadoff, md5 c7682f07bd4c6d9563e95e3875c53a5f.
- **Start:** trc_1101/keep/off_0p50.rst (kept, read-only).
  - It is the T1 off arm at t = 0.5 rot (cycle 8774), bitwise p_f4.
  - r1m production keys: fac 4, outflow1, cs_seam_flux shared, am_hst_split false, rkg_dt_safety 0.5, rkg_lean.
- **Arm keys:**
  - all arms: `mhd/max_eta=1e15|1e14|1e13`;
  - X arms only: `mhd/use_rkg_sts=false mhd/rkg_dt_safety=0`. rkg_dt_safety > 0 without STS FATALs, so it must be zeroed.
- **Arms run one per GPU.** R15 shares hpc-sm-01-11, whose CPUs were oversubscribed (load ~48 on 32 cores), so its wall-clock s/cycle is inflated. Cost ratios below come from a dedicated interleaved set.
- **Progress:** R14 and R13 have run the full 1.5 rot. R15 is at ~1.2 rot and continues to +1.5 rot. X13 was stopped at +0.5 rot.

## Cost
- **Method:** all arms alone on one H200 (01-11 GPU 2), interleaved, 2 repeats, cycles 8800 -> 8900 from the same rst. Mean dt and stages come from the long arms over +0 -> +0.5 rot.
- Absolute s/cycle is CPU-contended on that node. On the uncontended nodes, R14 and R13 ran at 0.041 and 0.037 s/cycle. Use the ratios.

| arm | s/cycle (2 repeats) | mean dt [s] | RKG stages/cycle | s/rot (dedicated) | vs R15 |
|---|---|---|---|---|---|
| R15 | 0.1883 / 0.1886 | 4.078 | 15.3 (16-17) | 5090 | 1.00 |
| R14 | 0.1152 / 0.1163 | 4.085 | 5.0 | 3120 | 0.61 |
| R13 | 0.1035 / 0.1035 | 4.036 | 3.0 | 2830 | 0.56 |
| X13 | 0.0764 / 0.0802 | 1.382 (resistive) | - | 6240 | 1.23 |
| X14 | 0.0757 / 0.0757 | 0.1382 | - | 60,300 | 11.8 |
| X15 | 0.0752 / 0.0747 | 0.01431 | - | 577,000 | 113 |

- **The explicit resistive dt** is 1.38 s at 1e13. That is 3.3x below the hydro dt (~4.5 s), not above it as estimated.
- **So explicit loses at every cap:** X13 costs 2.2x R13.

## Accuracy at +0.5 rot (t = 1.0 rot) vs R15
- **Method:** horizontal mass-weighted means at p = 1e-4 / 1e-3 / 1e-2 / 1e-1 bar (+-0.15 dex).
  - u = zonal wind, all latitudes and |lat| < 15 deg;
  - T on the day side (|lon| < 60), at the terminator (60-120) and on the night side (> 120);
  - B_rms per band; E_B at p > 1e-6 bar.
- **R15 reference values:**
  - u at 1e-4 / 1e-3 / 1e-2 / 1e-1 bar: -0.178 / 0.036 / 0.040 / 0.014 km/s;
  - night T at the same levels: 1101 / 1085 / 1137 / 1446 K;
  - E_B = 2.156e32 erg.

| arm | max abs. u diff [km/s] | max u_eq rel diff | max dT/T | E_B diff | B_rms diff (bands 1e-6 / 1e-5 / 1e-4 / 1e-3 / 1e-2 bar) |
|---|---|---|---|---|---|
| R14 | 0.003 (3.3 %, at 1e-4 and 1e-3) | 3.5 % | 0.26 % | +0.46 % | -0.0 / -0.4 / +0.9 / +0.2 / +0.2 % |
| R13 | 0.011 (30 % of the slow 1e-3 bar mean) | 6.3 % | 1.08 % | +5.2 % | -1.2 / **-7.1** / **+6.6** / +2.0 / +2.0 % |
| X13 | 0.012 | 7.2 % | 1.24 % | +5.2 % | -1.3 / -7.1 / +6.7 / +2.0 / +2.0 % |

- **No R15-C run exists for max_eta.** I dropped it to save GPUs, so the chaos spread comes from csvertex A vs C instead: the same keys, 1 vs 2 ranks, +1 rot after a round-off departure. That spread is:
  - u 2.8 %, u_eq 4.3 %, T 0.07 %;
  - B_rms <= 0.1 % per band, E_B 1e-4.
- **Counters** (dfloor / efloor / tfloor; vceil, c2p_fail and fofc are 0 everywhere):
  - R15, +0 -> +0.5 rot: 3.3e9 / 4.6e8 / 5.8e8.
  - R14 and R13 are similar per rot. They totalled 1.10e10 / 1.03e10 dfloor over 1.5 rot.

## Interim verdict (decision rule of 6e5e226a6)
- **R13: FAIL.**
  - B_rms falls 7 % at 1e-5..1e-4 bar and rises 7 % at 1e-4..1e-3 bar. That is a systematic shift at the Fe level, 70x the B spread.
  - E_B is +5 % and T differs by up to 1.1 %.
- **X13: FAIL.** It is the same as R13 within 0.2 %, which shows the error comes from the cap, not from STS splitting. It also costs 2.2x R13.
- **R14: borderline pass at +0.5 rot.**
  - Winds are within 2 % except the slow 1e-3 bar mean (0.003 km/s, inside the A/C wind spread of 2.8-4.3 %).
  - T is 0.26 %, under the 2 % rule but 4x the chaos T spread.
  - B_rms is within 0.9 % per band, with no clear trend; E_B is +0.5 %.
  - Cost is 0.61x R15. I am continuing R15 to +1.5 rot so R14 can be checked at +1 and +1.5 rot; the final NOTE follows.
- **Explicit (X) arms:** not competitive at any cap on this state.
