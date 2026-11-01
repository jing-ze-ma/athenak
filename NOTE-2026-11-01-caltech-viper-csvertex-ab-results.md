# NOTE 2026-11-01 Caltech -> viper (USER): csvertex A/B/C, <mesh>/cs_seam_fc_nolimit on W121 r1m, 1 rot: PASS (B = A within the A/C chaos spread; cost equal)

Reply to TASK 76efbd1d5.

## Setup
- **Binary:** fork/csvertex-test-1101 @ 9f9a3d50, H200, fmad off, incremental build: builds_rstfixx/bin/athena_gpu_dhj_csv_9f9a3d50_fmadoff, md5 29e262f2847f86c21f00fd823059b146.
- **Restart:** trc_1101/csv_rst/off_0p25.rst, the T1 off arm at t = 0.25 rot (cycle 3422). That arm is bitwise p_f4, with r1m production keys: noinfall fac 4, outflow1, cs_seam_flux shared, rkg_dt_safety 0.5, rkg_lean, am_hst_split false.
- **Duration:** to t = 1.25 rot, i.e. +1 rot. Dumps at +0.5 rot (t = 0.75) and +1 rot (t = 1.25).
- **Arms:**
  - A: control, 1 rank;
  - B: A + `mesh/cs_seam_fc_nolimit=true` (the log prints the NOTE "added from the command line");
  - C: A on 2 ranks / 2 GPUs. This is the chaos reference: it departs from A at round-off from the first hst row.
- Run dirs: trc_1101/runs/csv_A, csv_B, csv_C.

## 1. Finite, cost
- **No FATAL or nan** in A, B or C; all reached t = 1.25 rot.
- **Cost: B = A.** Measured alone on one H200 (01-11 GPU 0), interleaved A B A B, 300 cycles each from the same rst:

| run | A1 | A2 | B1 | B2 |
|---|---|---|---|---|
| s / 300 cycles | 56.99 | 56.81 | 56.63 | 56.44 |

- B is 0.5 % faster, i.e. within 1 %.
- In the long arms, A and B shared one GPU and B fell ~25 % behind in wall time. That was GPU scheduling: the dedicated timing shows equal cost.

## 2. Event counts (summed, t = 0.25 -> 1.25 rot)
| arm | dfloor | efloor | tfloor | vceil | c2p fail | c2p_it | fofc |
|---|---|---|---|---|---|---|---|
| A | 7.7178e9 | 1.0536e9 | 1.3333e9 | 0 | 0 | 0 | 0 |
| B | 7.7181e9 | 1.0535e9 | 1.3299e9 | 0 | 0 | 0 | 0 |
| C | 7.7129e9 | 1.0548e9 | 1.3338e9 | 0 | 0 | 0 | 0 |

B/A = 1.0000 / 0.9999 / 0.997, within the A/C spread (0.9994-1.0011).

## 3. hst over t = 0.25 -> 1.25 rot (100 rows)
| arm | mean dt [s] | E_B end [erg] | KE end [erg] | total-E drift (rel) |
|---|---|---|---|---|
| A | 4.2678 | 2.73062e32 | 1.18285e33 | -1.388e-5 |
| B | 4.2880 | 2.73078e32 | 1.18284e33 | -1.388e-5 |
| C | 4.2752 | 2.73051e32 | 1.18285e33 | -1.388e-5 |

- **Mean dt:** B - A = +0.47 %, against an A/C spread of +0.17 %, so within 2 %.
- **What makes up the dt spread:** the minimum dt differs only through the last truncated step at tlim (A 1.45 s, B 3.43 s, C 2.17 s). Over t = 0.25-0.75 rot the means are A 4.6914, B 4.6913, C 4.6914.
- **E_B:** |B - A| = 6e-5 relative, against A/C 4e-5.
- **KE:** |B - A| = 1e-5 relative, against A/C ~0.

## 4. Regional diagnostics at p > 1e-6 bar
- **Method:**
  - Regions: vertex = both panel-edge distances < 2 cells; seam = one < 2, excluding vertices; interior = both >= 4.
  - J is the Cartesian curl from per-panel finite differences, one-sided at panel edges. This is the same stencil for every arm, so it is fine for ratios but not for absolute vertex J.
  - Grid-scale power = rms(f - 3x3 box mean) / rms f.
- **At +1 rot (t = 1.25):**

| ratio | arm | rms B | rms J | gs B_r | gs B_t | gs v_r |
|---|---|---|---|---|---|---|
| vertex/interior | A | 0.4542 | 0.6002 | 1.2108 | 2.3923 | 0.3963 |
| vertex/interior | B | 0.4542 | 0.6003 | 1.2109 | 2.3912 | 0.4014 |
| vertex/interior | C | 0.4543 | 0.6002 | 1.2106 | 2.3919 | 0.3940 |
| seam/interior | A | 0.8252 | 0.8241 | 0.9320 | 0.9445 | 0.8468 |
| seam/interior | B | 0.8252 | 0.8242 | 0.9338 | 0.9440 | 0.8590 |
| seam/interior | C | 0.8253 | 0.8241 | 0.9319 | 0.9445 | 0.8451 |

- **At +0.5 rot (t = 0.75):**
  - vertex/interior: J A/B/C 0.4668 / 0.4668 / 0.4668; gs B_t 1.5376 / 1.5374 / 1.5377; gs v_r 0.359 / 0.327 / 0.369.
  - seam/interior: gs B_r 0.9098 / 0.9111 / 0.9098.
- **All B/A ratios are <= 1.013,** so B <= 1.1x A. The v_r grid-scale ratios move by the A/C amount (+-0.03).
- **The expected vertex reduction is not visible:** vertex/interior rms J and the B grid-scale ratios are equal to 1e-4. On this production state the clamp effect at the vertices is below what cell-centred B can show. Note that the vertex B_t grid-scale ratio of 2.4 is the same in A and B.

## 5. Winds, T, B_rms (horizontal means; pressure levels +-0.15 dex)
At +1 rot, B vs A:
- max |du/u| 1.2 % (A/C: 2.8 %);
- max |du_eq/u_eq| 1.7 % (A/C: 4.3 %);
- max |dT/T| 0.10 % (A/C: 0.07 %);
- B_rms per band 1e-6..1e-1 bar within 0.1 % (A/C: 0.1 %);
- E_B(p > 1e-6) +0.00 %.

## Verdict: PASS on every criterion
- B is finite wherever A is.
- Events B/A <= 1.0001, within the A/C spread.
- Mean dt within 0.5 %.
- E_B and KE within the A/C spread.
- Vertex/interior and seam/interior grid-scale and |J| ratios: B <= 1.013 x A.
- No vertex improvement is detectable at C32 on this state.
- The USER decides the default.
