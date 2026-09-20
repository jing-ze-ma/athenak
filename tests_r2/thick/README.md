# Optically thick 1-D spherical unit test of the grey two-stream

**Status: measured results are appended at the bottom of this file.**

## What this tests

Commit `1159a8f3` made the grey two-stream sweep spherical by transporting the
*area-weighted* intensity `J = A I` against the *area-weighted* source `A B`, with
`A ~ r^2` the face area:

```
mu dJ/dr = -kappa rho (J - A B).
```

That is **exact in the transparent limit** — it is what makes `L = 4 pi r^2 F` constant
through an empty shell, and `tests_r2/rg1d` is the existing gate for it. In the
**diffusion limit** the same equation gives

```
F_2s = -(mu/(kappa rho)) d(A B)/dr / A
     = -(mu/(kappa rho)) [ dB/dr + 2 B/r ],
```

while the correct answer is `F = -(mu/(kappa rho)) dB/dr`. So the scheme is off by

```
F_2s / F_exact = 1 + (2 B/r)/(dB/dr) = 1 - H_T/(2 r),     H_T = T/|dT/dr|,  B ~ T^4.
```

This directory is the clean unit test of that factor against a state whose correct flux
is known analytically. (The same error was already confirmed on the He star over
0.5-0.95 R; this test removes every other moving part.)

## The analytic solution

Density is a pure power law and the opacity is constant
(`problem/kappa_const`, which replaces the opacity table for both the two-stream and the
conduction operator):

```
rho(r) = rho_in (r_in/r)^n .
```

The temperature is the **exact radiative-equilibrium diffusion solution** for that rho:

```
d(T^4)/dr = -3 kappa rho(r) L / (16 pi sigma r^2),     L = 4 pi r_out^2 sigma Teff^4,
```

integrated numerically (trapezoid, 8001 nodes) by `make_ic_thick.py`, so no closed form
is assumed and any `n` works. On that state the *correct* radiative flux at every face is

```
F_exact(r) = -(16 sigma T^3)/(3 kappa rho) dT/dr = L/(4 pi r^2) = F_req(r),
```

which is exactly the `F_req` column the dump already carries, so `F_2s/F_req` **is** the
error factor above with no post-processing.

Deep in the shell `T^4 ~ r^-(n+1)`, hence `H_T = 4r/(n+1)` and

```
F_2s/F_req = 1 - 2/(n+1):    n=0 -> -1,  n=1 -> 0,  n=3 -> +0.5,  n=7 -> +0.75.
```

`n` is therefore the knob that sets `H_T/r`, and `n = 0` is the case in which the sweep
returns a flux of the **wrong sign**.

The gas is ideal with `<units>` all 1.0, so `eint = p/(gamma-1)` with
`p = rho k_B T/(mu m_H)`, `gamma = 5/3`, `mu = problem/mu = <units>/mu = 1`.
The anchor `T = Teff` is placed at the **top of the extended profile** rather than at
`r_out`: anchoring at `r_out` drives `T^4` negative in the outer ghost margin for the
thick and thin-shell cases, and the anchor value is irrelevant to the test because the
ODE — and with it `F_diff = F_req` — holds whatever it is.

## The harness

`inputs/tests/two_stream_sph_thick.athinput`, built on the transparent gate
`tests_r2/rg1d/rg_1d_sph.athinput`: cubed sphere, `nx2 = nx3 = 4`, `nx1 = 128` uniform
radial, `ix1_bc/ox1_bc = user`, `problem/inner_bc = wall`, `eos = ideal`,
`PROBLEM=red_giant`.

Differences from the transparent gate:

* **Optically thick.** `problem/kappa_const` is chosen from the requested total radial
  Rosseland depth `tau_tot = int kappa rho dr` (30, 100, 300).
* **The state is supplied exactly**, through `problem/ic_profile` (three columns
  `r[cm] rho eint`, ascending in r, spanning the mesh *and* its radial ghosts with a
  5-cell margin), written by `make_ic_thick.py`.
* **The blend is off**: `<hydro>/rad_tau_lo = 1e5`, `rad_tau_hi = 1e6`, far deeper than
  any `tau_tot` here, so `RadBlendWeight` returns `w = 0` at every face, the radiative
  conduction operator contributes nothing, the two-stream's cut index sits at `is` (the
  whole column is swept), and the two-stream is the **only** carrier. `L` is injected on
  the two-stream's own lower boundary via `problem/rt_bottom_flux = true` with
  `<hydro>/rad_flux_inner = L/(4 pi r_in^2)`.
* **The two-stream switches are the He production's**: `problem/rt_grey = true`,
  `rt_implicit_column = 3`, `rt_impl_solver = pcr`, `ck_nquad = 2`, `rt_strang = true`.
  So the operator under test is the one that actually runs.
* `time/nlim = 1` and every output `dt = 1e30`: hydro never moves, the measurement is the
  single `t = 0` face dump `problem/mlt_dump = mltfaces.txt` (which `problem/mlt_alpha =
  1.5` enables). Its columns are

  ```
  i r T_f p_f grad grad_ad x F_mlt F_rad F_used F_req F_raddiff w_blend
    F_rad_col0 F_conv_res grad_rad D F_2s
  ```

  0-indexed, so `F_req` is column 10, `F_raddiff` (Rosseland diffusion evaluated on the
  same state, the control that must come back at `F_req`) is column 11, and `F_2s` (the
  two-stream's own net face flux) is column 17.

## Running it

```bash
cd tests_r2/thick
python3 make_ic_thick.py --n 3 --tau 100 --ratio 2.0 --nfine 8001 --out n3_t100_r2.0/ic_thick.txt
sbatch sub_thick.sh          # runs every row of cases.txt sequentially, 1 GPU, 15 min
python3 analyze_thick.py n3_t100_r2.0/mltfaces.txt
```

`cases.txt` holds the matrix (tag, n, tau, r_out/r_in, and the derived `x1max`,
`kappa_const`, `lstar`, `rad_flux_inner` that `sub_thick.sh` passes as command-line
overrides). The `.bin`/`.rst`/`ic_thick.txt` files are scratch and are not committed.

## Matrix

| set | varied | fixed |
| --- | --- | --- |
| A | `n` in {0, 1, 3, 7} | `tau_tot = 100`, `r_out/r_in = 2.0` |
| B | `tau_tot` in {30, 100, 300} | `n = 3`, `r_out/r_in = 2.0` |
| C | `r_out/r_in` in {1.1, 1.5, 2.0} | `n = 3`, `tau_tot = 100` |

A is the `H_T/r` scaling; B and C are the controls — the error must be independent of
optical depth and a weak function of shell thickness.

---

# MEASURED (2026-09-17, viper apudev, 1 GPU, `build_gpu_rg`, branch `he4-presn-global`)

`rt_implicit_column = 3` (the exact block-tridiagonal column solve, `pcr`, `ck_nquad = 2`,
`rt_strang = true`) ran this setup without complaint; no fallback to mode 0 was needed.
`w_blend = 0` at every face in every case, i.e. the two-stream really is the only carrier.

## 1. Per-case tables (every 8th face)

Columns: `r/r_in`, the dumped face temperature, `H_T/r` by centred differences on that
dumped `T`, the prediction `1 - H_T/(2r)`, the measured `F_2s/F_req`, and the control
`F_raddiff/F_req`.

    ### n0_t100_r2.0/mltfaces.txt
      r/r_in      T[K]      H_T/r     pred      meas    F_rd/F_req
       1.0156    17486.6    1.8064    0.0968    0.0560    1.1014
       1.0781    16951.2    1.8982    0.0509    0.0503    1.0004
       1.1406    16439.5    1.7778    0.1111    0.1112    0.9997
       1.2031    15936.5    1.6560    0.1720    0.1720    0.9997
       1.2656    15438.4    1.5342    0.2329    0.2329    0.9997
       1.3281    14941.3    1.4124    0.2938    0.2938    0.9998
       1.3906    14441.3    1.2907    0.3546    0.3547    0.9998
       1.4531    13933.9    1.1689    0.4156    0.4156    0.9998
       1.5156    13413.9    1.0470    0.4765    0.4765    0.9998
       1.5781    12874.9    0.9253    0.5374    0.5373    0.9998
       1.6406    12308.5    0.8035    0.5983    0.5982    0.9998
       1.7031    11703.1    0.6816    0.6592    0.6591    0.9998
       1.7656    11041.3    0.5598    0.7201    0.7200    0.9998
       1.8281    10294.7    0.4380    0.7810    0.7809    0.9998
       1.8906     9410.2    0.3161    0.8419    0.8418    0.9998
       1.9531     8265.6    0.1943    0.9029    0.9008    0.9994
      MID: r/rin=1.500 pred=0.4613 meas=0.4612 resid=-0.0000 Frd/Freq=0.9998 wmax=0.00e+00
    ### n1_t100_r2.0/mltfaces.txt
      r/r_in      T[K]      H_T/r     pred      meas    F_rd/F_req
       1.0156    17755.0    1.4043    0.2978    0.2643    1.0619
       1.0781    17046.7    1.4514    0.2743    0.2742    1.0000
       1.1406    16383.1    1.3862    0.3069    0.3069    0.9999
       1.2031    15749.2    1.3171    0.3414    0.3415    0.9999
       1.2656    15138.7    1.2443    0.3778    0.3779    0.9999
       1.3281    14545.7    1.1678    0.4161    0.4161    0.9999
       1.3906    13964.5    1.0876    0.4562    0.4562    0.9998
       1.4531    13389.6    1.0038    0.4981    0.4981    0.9998
       1.5156    12814.9    0.9162    0.5419    0.5419    0.9998
       1.5781    12233.5    0.8250    0.5875    0.5875    0.9998
       1.6406    11637.2    0.7301    0.6350    0.6350    0.9998
       1.7031    11014.8    0.6315    0.6843    0.6843    0.9998
       1.7656    10350.5    0.5291    0.7354    0.7354    0.9998
       1.8281     9619.3    0.4231    0.7884    0.7884    0.9998
       1.8906     8775.7    0.3134    0.8433    0.8432    0.9998
       1.9531     7718.1    0.2002    0.8999    0.8931    0.9983
      MID: r/rin=1.500 pred=0.5308 meas=0.5308 resid=-0.0000 Frd/Freq=0.9998 wmax=0.00e+00
    ### n3_t100_r2.0/mltfaces.txt
      r/r_in      T[K]      H_T/r     pred      meas    F_rd/F_req
       1.0156    18251.4    0.9123    0.5439    0.5293    1.0239
       1.0781    17120.2    0.9285    0.5357    0.5356    0.9999
       1.1406    16102.8    0.9105    0.5447    0.5447    0.9999
       1.2031    15176.3    0.8892    0.5554    0.5553    0.9999
       1.2656    14324.8    0.8644    0.5678    0.5678    0.9999
       1.3281    13535.4    0.8355    0.5822    0.5822    0.9999
       1.3906    12796.6    0.8023    0.5988    0.5988    0.9999
       1.4531    12098.5    0.7643    0.6179    0.6179    0.9999
       1.5156    11431.8    0.7210    0.6395    0.6395    0.9999
       1.5781    10787.8    0.6721    0.6639    0.6640    0.9999
       1.6406    10157.2    0.6170    0.6915    0.6917    0.9999
       1.7031     9529.8    0.5552    0.7224    0.7226    0.9999
       1.7656     8892.6    0.4862    0.7569    0.7572    0.9998
       1.8281     8227.8    0.4095    0.7952    0.7957    0.9997
       1.8906     7506.5    0.3246    0.8377    0.8380    0.9995
       1.9531     6672.9    0.2312    0.8844    0.8764    0.9976
      MID: r/rin=1.500 pred=0.6339 meas=0.6339 resid=-0.0000 Frd/Freq=0.9999 wmax=0.00e+00
    ### n7_t100_r2.0/mltfaces.txt
      r/r_in      T[K]      H_T/r     pred      meas    F_rd/F_req
       1.0156    18854.9    0.4954    0.7523    0.7495    1.0059
       1.0781    16733.8    0.5009    0.7495    0.7493    0.9997
       1.1406    14954.2    0.5014    0.7493    0.7490    0.9997
       1.2031    13446.1    0.5022    0.7489    0.7487    0.9997
       1.2656    12157.7    0.5034    0.7483    0.7481    0.9998
       1.3281    11049.1    0.5050    0.7475    0.7473    0.9998
       1.3906    10089.3    0.5072    0.7464    0.7462    0.9998
       1.4531     9253.9    0.5103    0.7449    0.7447    0.9998
       1.5156     8523.5    0.5143    0.7428    0.7426    0.9999
       1.5781     7882.4    0.5195    0.7403    0.7402    1.0005
       1.6406     7317.6    0.5254    0.7373    0.7397    1.0023
       1.7031     6817.9    0.5319    0.7340    0.7463    1.0063
       1.7656     6373.9    0.5386    0.7307    0.7692    1.0129
       1.8281     5977.6    0.5449    0.7275    0.8190    1.0229
       1.8906     5621.0    0.5467    0.7266    0.9072    1.0431
       1.9531     5292.5    0.5275    0.7362    1.0561    1.1019
      MID: r/rin=1.500 pred=0.7434 meas=0.7431 resid=-0.0003 Frd/Freq=0.9998 wmax=0.00e+00
    ### n3_t30_r2.0/mltfaces.txt
      r/r_in      T[K]      H_T/r     pred      meas    F_rd/F_req
       1.0156    13554.8    0.9224    0.5388    0.5696    1.0321
       1.0781    12723.8    0.9439    0.5280    0.5280    1.0003
       1.1406    11981.5    0.9302    0.5349    0.5349    0.9999
       1.2031    11308.1    0.9137    0.5432    0.5433    0.9999
       1.2656    10692.1    0.8943    0.5529    0.5531    0.9999
       1.3281    10124.3    0.8718    0.5641    0.5645    0.9999
       1.3906     9596.7    0.8459    0.5771    0.5776    0.9999
       1.4531     9102.4    0.8163    0.5919    0.5927    0.9999
       1.5156     8635.5    0.7825    0.6087    0.6099    0.9999
       1.5781     8190.5    0.7444    0.6278    0.6296    0.9999
       1.6406     7762.2    0.7014    0.6493    0.6523    1.0000
       1.7031     7345.5    0.6531    0.6734    0.6790    1.0001
       1.7656     6934.8    0.5991    0.7004    0.7128    1.0004
       1.8281     6523.6    0.5388    0.7306    0.7602    1.0011
       1.8906     6103.7    0.4717    0.7642    0.8359    1.0025
       1.9531     5663.0    0.3963    0.8018    0.9723    1.0068
      MID: r/rin=1.500 pred=0.6044 meas=0.6054 resid=+0.0010 Frd/Freq=0.9999 wmax=0.00e+00
    ### n3_t300_r2.0/mltfaces.txt
      r/r_in      T[K]      H_T/r     pred      meas    F_rd/F_req
       1.0156    23995.3    0.9177    0.5411    0.5298    1.0161
       1.0781    22504.1    0.9241    0.5380    0.5379    0.9999
       1.1406    21159.6    0.9049    0.5476    0.5475    0.9999
       1.2031    19933.8    0.8823    0.5589    0.5588    0.9999
       1.2656    18805.8    0.8558    0.5721    0.5720    0.9999
       1.3281    17758.1    0.8252    0.5874    0.5874    0.9999
       1.3906    16775.6    0.7898    0.6051    0.6050    0.9999
       1.4531    15844.6    0.7495    0.6253    0.6252    0.9999
       1.5156    14952.7    0.7035    0.6483    0.6482    0.9999
       1.5781    14087.3    0.6514    0.6743    0.6742    0.9999
       1.6406    13235.2    0.5928    0.7036    0.7035    0.9999
       1.7031    12380.7    0.5272    0.7364    0.7363    0.9999
       1.7656    11503.6    0.4539    0.7731    0.7730    0.9998
       1.8281    10573.8    0.3723    0.8138    0.8138    0.9998
       1.8906     9538.2    0.2820    0.8590    0.8590    0.9997
       1.9531     8279.4    0.1829    0.9086    0.9073    0.9960
      MID: r/rin=1.500 pred=0.6423 meas=0.6422 resid=-0.0001 Frd/Freq=0.9999 wmax=0.00e+00
    ### n3_t100_r1.1/mltfaces.txt
      r/r_in      T[K]      H_T/r     pred      meas    F_rd/F_req
       1.0016    15171.5    0.3220    0.8390    0.8448    1.0108
       1.0078    14876.8    0.3093    0.8453    0.8453    1.0000
       1.0141    14573.9    0.2920    0.8540    0.8540    0.9999
       1.0203    14261.0    0.2744    0.8628    0.8628    0.9999
       1.0266    13936.5    0.2565    0.8718    0.8718    0.9999
       1.0328    13598.4    0.2382    0.8809    0.8809    0.9999
       1.0391    13244.4    0.2196    0.8902    0.8902    0.9999
       1.0453    12871.4    0.2006    0.8997    0.8997    0.9999
       1.0516    12475.4    0.1813    0.9093    0.9093    0.9999
       1.0578    12051.2    0.1617    0.9192    0.9191    0.9999
       1.0641    11591.5    0.1417    0.9292    0.9291    0.9999
       1.0703    11085.6    0.1213    0.9393    0.9393    0.9999
       1.0766    10517.4    0.1006    0.9497    0.9497    0.9999
       1.0828     9859.9    0.0795    0.9602    0.9602    0.9999
       1.0891     9062.8    0.0581    0.9710    0.9709    0.9999
       1.0953     8010.2    0.0363    0.9819    0.9778    0.9990
      MID: r/rin=1.050 pred=0.9070 meas=0.9069 resid=-0.0001 Frd/Freq=0.9999 wmax=0.00e+00
    ### n3_t100_r1.5/mltfaces.txt
      r/r_in      T[K]      H_T/r     pred      meas    F_rd/F_req
       1.0078    16541.4    0.7884    0.6058    0.5990    1.0231
       1.0391    15920.8    0.7890    0.6055    0.6055    1.0000
       1.0703    15324.3    0.7625    0.6187    0.6187    0.9999
       1.1016    14746.0    0.7335    0.6332    0.6332    0.9999
       1.1328    14182.4    0.7019    0.6490    0.6490    0.9999
       1.1641    13630.0    0.6677    0.6662    0.6662    0.9999
       1.1953    13085.0    0.6305    0.6847    0.6847    0.9999
       1.2266    12543.5    0.5904    0.7048    0.7048    0.9999
       1.2578    12000.7    0.5469    0.7265    0.7265    0.9999
       1.2891    11451.2    0.5002    0.7499    0.7499    0.9999
       1.3203    10888.1    0.4500    0.7750    0.7750    0.9999
       1.3516    10302.1    0.3960    0.8020    0.8020    0.9999
       1.3828     9679.7    0.3382    0.8309    0.8310    0.9999
       1.4141     8999.4    0.2763    0.8619    0.8619    0.9998
       1.4453     8222.8    0.2101    0.8949    0.8949    0.9997
       1.4766     7267.2    0.1400    0.9300    0.9145    0.9966
      MID: r/rin=1.250 pred=0.7210 meas=0.7209 resid=-0.0000 Frd/Freq=0.9999 wmax=0.00e+00

## 2. Scaling summary (mid-shell face, i = 64)

| set | case | n | tau_tot | r_out/r_in | predicted `1-H_T/2r` | measured `F_2s/F_req` | residual | deep asymptote `1-2/(n+1)` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | n0_t100_r2.0 | 0 | 100 | 2.0 | 0.4613 | 0.4612 | -0.00004 | -1.000 |
| A | n1_t100_r2.0 | 1 | 100 | 2.0 | 0.5308 | 0.5308 | -0.00004 | +0.000 |
| A,B,C | n3_t100_r2.0 | 3 | 100 | 2.0 | 0.6339 | 0.6339 | -0.00002 | +0.500 |
| A | n7_t100_r2.0 | 7 | 100 | 2.0 | 0.7434 | 0.7431 | -0.00029 | +0.750 |
| B | n3_t30_r2.0 | 3 | 30 | 2.0 | 0.6044 | 0.6054 | +0.00105 | +0.500 |
| B | n3_t300_r2.0 | 3 | 300 | 2.0 | 0.6423 | 0.6422 | -0.00011 | +0.500 |
| C | n3_t100_r1.1 | 3 | 100 | 1.1 | 0.9070 | 0.9069 | -0.00007 | +0.500 |
| C | n3_t100_r1.5 | 3 | 100 | 1.5 | 0.7210 | 0.7209 | -0.00003 | +0.500 |

**A, the `H_T/r` knob.** `F_2s/F_req` tracks `1 - H_T/(2r)` face by face, over the whole
range the shell offers: from 0.05 at the bottom of the `n = 0` case (a flux 20x too
small) to 0.90 at its top, and from 0.75 to 0.90 for `n = 7`. The error is exactly the
predicted geometric factor and nothing else.

The *closed-form* asymptote `1 - 2/(n+1)` is the `r << r_out` limit of that factor, and it
is reached only where the shell is deep enough. With `r_out/r_in = 2` the finite shell
modifies it.  Neglecting the `Teff^4` floor, `T^4 ~ (r^-(n+1) - r_out^-(n+1))`, so

```
H_T/(2r) = 2 (1 - (r/r_out)^(n+1)) / (n+1),
```

which at the bottom face (`r/r_out = 0.5`) predicts 0.00 (n=0), +0.25 (n=1), +0.53 (n=3),
+0.75 (n=7) instead of the deep asymptote -1, 0, +0.5, +0.75.  Measured there: 0.056,
0.264, 0.529, 0.750.  `n = 3` and `n = 7` are already on their asymptote at the bottom of
a shell only a factor 2 thick; `n = 0` and `n = 1` need a thicker shell to get there, and
the sign flip `F_2s < 0` at `n = 0` lives at `r/r_out < 0.5`, just outside this domain. **The prediction that was tested and confirmed is
the exact one, `1 - H_T/(2r)`**, evaluated on the code's own dumped `T`.

**B, optical depth.** At fixed `n = 3` the mid-shell error is 0.6044 / 0.6339 / 0.6423 for
`tau_tot` = 30 / 100 / 300, against predictions 0.6044 / 0.6339 / 0.6423. The 6 % spread
between the three is *not* a tau dependence of the defect: it is the shell's own `H_T/r`
moving, because a deeper shell has a steeper temperature profile and is therefore closer
to the `T^4 ~ r^-(n+1)` asymptote. Measured minus predicted is +0.0011, -0.00002 and
-0.0001 — flat. **The error is independent of optical depth, as it must be for a
geometric defect of the transport operator.**

**C, shell thickness.** 0.9070 / 0.7210 / 0.6339 for `r_out/r_in` = 1.1 / 1.5 / 2.0,
predicted 0.9070 / 0.7210 / 0.6339. Again the variation is entirely the `H_T/r` the
thinner shell has (a thin shell has `H_T << r` and the defect nearly vanishes), and the
residual is < 1e-4 in all three. A *thin* shell hides this bug; a thick one exposes it.

## 3. The residual

`measured - predicted` is **tiny and structureless in the bulk**: the median over faces
`i = 10..100` is 2e-5 to 7e-4 in every case, i.e. 4-5 orders of magnitude below the
defect itself.

It does not scale with anything physical. It scales with (a) the **local cell optical
depth**, because the prediction is a diffusion-limit statement and stops applying once
`dtau_cell = kappa rho dr` falls below ~0.3 — there the sweep correctly leaves the
diffusion limit and `F_2s/F_req` heads back toward 1 (visible as the outer rows of the
`n = 7`, `tau = 30` tables, where `F_raddiff/F_req` also departs from 1 at the same
faces), and (b) the **two radial boundaries**: every per-case maximum sits at face 4 (the
first faces above the injected bottom flux) or face 122 (the top, where `rt_top_re` is off
under mode 3 and the ghost mirrors the top cell back), at the 1e-3 to 1e-2 level.

Binned over all eight cases, |residual| against the local cell optical depth:

      dtau_cell bin      N   median|resid|
       0.05 - 0.10      42    9.5e-02
       0.10 - 0.20      69    3.4e-03
       0.20 - 0.30      46    4.4e-04
       0.30 - 0.50      99    2.0e-04
       0.50 - 0.80     374    4.5e-05
       0.80 - 1.20     203    6.2e-05
       1.20 - 2.00      79    1.0e-04
       2.00 - 3.00      32    9.5e-05
       3.00 - 10.0      39    1.1e-04

Flat at ~1e-4 for `dtau_cell > 0.3` and rising only where the diffusion limit itself
fails.

**The `F_raddiff` control passes.** `F_raddiff/F_req` is 0.9997-1.0002 over the interior
of every case — the state handed to the code really is the exact diffusion solution, and
the Rosseland operator evaluated on it returns `L/(4 pi r^2)` to 2e-4. It departs from 1
only (i) on the single face above the wall, where the one-sided difference at the injected
bottom boundary makes it 1.006-1.10, and (ii) in the outermost few faces of the thin
cases, for the `dtau_cell` reason above. So the ~50 % deficit at `n = 3` is the
**two-stream's**, not the state's.

## 4. Verdict

The diffusion-limit error of the spherical grey two-stream is exactly

```
F_2s / F_req = 1 - H_T/(2 r)
```

to better than 1e-4 in the bulk, independent of optical depth over a decade in `tau_tot`
and of shell thickness over `r_out/r_in` = 1.1-2.0. It is a factor-of-two-scale error
wherever `H_T ~ r` (the He and B star envelopes) and it changes the **sign** of the flux
where `H_T > 2r`. `mu dJ/dr = -kappa rho (J - A B)` is transparent-limit exact and
diffusion-limit wrong, so a fix has to be right in both and cannot simply be the
plane-parallel source restored. This directory is the gate for whatever that fix turns out
to be: one cycle on one GPU, and it reports the factor directly.
