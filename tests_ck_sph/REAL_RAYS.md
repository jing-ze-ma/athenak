# The stellar beam near the terminator, checked against real rays

**Question.** Near the terminator, a grazing stellar ray passes through neighbouring columns, not only
the target column. Neither beam rule models that: the OLD rule (`ck_beam_sph = false`) and the NEW
rule (`ck_beam_sph = true`, used in prod4) both fill every shell with the target column's own opacity.
The NEW rule corrects the geometric error only. This note measures how the two errors compare
on a developed 3-D atmosphere.

**Answer.**
- **Dayside, 0 < mu0 < 0.45:** the NEW rule is closer to the real rays in every bin except the
  last one (0.40 to 0.45), where the two are equal to within 0.6 %. NEW over-deposits by 3 to 13 %.
  OLD under-deposits by 2 to 23 %, except in the first bin (0 to 0.05), where it over-deposits by 18 %.
- **Twilight ring, mu0 < 0:** OLD deposits nothing, which is wrong. NEW deposits beam power there, but
  **too little**: 0.99 of the real value near mu0 = 0, falling to 0.5 for mu0 between -0.25 and -0.45.
  In that ring the horizontal-inhomogeneity error (a factor of about 2) is larger than the geometry
  error NEW fixes.
- **Global budget:** the band -0.55 <= mu0 < 0.45 carries 25 % of the planet's beam power. In that
  band NEW is +1.7 % and OLD -16 % off the real rays. Planet total: NEW +0.4 %, OLD -4.8 %.

Figure: `REAL_RAYS.png`. Scripts: `rr_lib.py`, `rr_validate.py`, `real_rays.py`, `real_rays_report.py`.
Data: `real_rays.npz`. All paths are in `tests_ck_sph/`.

## Method

* **Dump.** `bench/cs_mhd_prod3/bin/dhj.mhd_w_bcc.00142.bin`, t = 283.28 rotations. The directory was
  only read. It is a cubed-sphere run with 6 x 32 x 32 x 128 cells. The stretched radial faces come
  from the input's `r_poly` map.
  - T and p are from the general-EOS inversion (`docs/handover/scripts/dhjcs.py`), using the table
    `bench/cs_ens/analysis/eos_table.txt` (X = 0.7381, Y = 0.2485, H2 + ionisation + metal ionisation).
  - That table's header does not record metal condensation, which prod3 has on. This is a small
    caveat on T.
* **Opacity.** `rr_lib.py` transcribes `src/utils/correlated_k.hpp` line for line, with all 88 chains
  (11 bands x 8 g-points):
  - the line table, bilinear in log kappa;
  - the continuum: FastChem composition, 4 CIA pairs, 4 Rayleigh species, H- bound-free and
    free-free;
  - the 6000 K blackbody band fractions;
  - the kernel's two deposit formulas.
  Because all 88 chains are used, the chains cover 100 % of the stellar flux.
* **Cut.** prod3 runs the tau blend, so `icut` is the first face whose Rosseland depth is below
  `rad_tau_hi` = 300. That depth is rebuilt from `ck_build_rosseland_table`, transcribed. Result:
  icut = 24 to 25 in every column, median r_cut = 1.089e10 cm. The 10 bar level would give 16.
* **Gate 0: the transcription against the kernel.** `rr_validate.py` recomputes Q_sw of the kernel's
  own column dumps `bd_{p979,p435,p102}_{old,new}` and `bd_{n102,n195,n513}_new`, with rho recovered
  from the dumped (p, T).
  - Column power, Python / code: 0.99998 to 1.00000.
  - Per cell: 0.9990 to 1.0001.
* **Rays.** The star lies along -x, which matches `CSCellAngles`: mu0 = -x-hat. For each face of a
  target column, the ray starts at r_face n-hat and runs straight toward the star.
  - **Segments.** The segment edges are every exact shell crossing plus a uniform 2e7 cm grid. That grid
    is 6 % or less of an angular cell, whose smallest size is about 3.5e8 cm.
  - **Opacity along the ray.** Each segment uses kappa rho of the cell that contains its midpoint
    (nearest cell in 3-D).
  - **Termination.** The ray stops at r_top. It is dark if it enters any cell below that column's own
    icut.
  - **Deposit.** Qb is the kernel's `ck_beam_sph` deposit (README_beam.md section 2), applied to these
    optical depths.
  - **Old rule.** OLD is the kernel's facsw path, including the ghost column above the top. That ghost
    column is approximated with the top cell's kappa p/g, because the dump has no ghost cells.
  - **Column power.** P = sum over cells of Qb V, in units of (1-A) F\* a_p^2 per steradian. Albedo
    cancels in every ratio.
* **Sample.** 160 columns: 20 mu0 bins of width 0.05 from -0.55 to +0.45. Each bin has 4 columns with
  lon > 0 and 4 with lon < 0, spread in |lat| from the equator to the highest available. OLD and NEW
  are also evaluated in all 6144 columns for the planet totals.
* **Gate 1: geometry control.** The same real-ray code, run through a horizontally uniform atmosphere
  built from the target column, reproduces NEW:
  - To **2.4e-9 per cell** and **7e-10 in column power**, in every cell except the top cell (i = 127).
  - For mu0 < 0 that top cell differs by 11 to 18 %, which is at most 1.4 % of column power. The
    difference is a kernel convention: the kernel sets tau(r_top) = 0, but for mu0 < 0 the true chord
    from r_top dips back into the atmosphere.
  - The real rays keep the honest top-face depth. The convention is part of NEW's error. At most it is
    7.7 % of NEW's column power, for mu0 = -0.535, where the column power is 1e-3 of the dayside value.

## Numbers (`real_rays_report.py`)

Column power ratios are per column; each bin has n = 8. The two "cell" columns give per-cell ratios
over cells with Qb_real > 1 % of that column's peak.

| mu0 bin | P_real | P_new/P_real med [min, max] | P_old/P_real med [min, max] | cell Qnew/Qreal med [p10, p90] | cell Qold/Qreal med [p10, p90] |
|---|---|---|---|---|---|
| -0.55 .. -0.50 | 3.9e-04 | 1.239 [0.427, 1.639] | 0 | 1.497 [0.170, 2.205] | 0 |
| -0.50 .. -0.45 | 8.3e-04 | 0.840 [0.397, 1.309] | 0 | 1.277 [0.209, 2.843] | 0 |
| -0.45 .. -0.40 | 5.2e-03 | 0.505 [0.421, 0.940] | 0 | 0.936 [0.216, 1.910] | 0 |
| -0.40 .. -0.35 | 1.05e-02 | 0.515 [0.456, 0.586] | 0 | 0.885 [0.237, 1.335] | 0 |
| -0.35 .. -0.30 | 2.14e-02 | 0.517 [0.488, 0.593] | 0 | 0.883 [0.280, 1.289] | 0 |
| -0.30 .. -0.25 | 3.44e-02 | 0.567 [0.526, 0.634] | 0 | 0.941 [0.336, 1.343] | 0 |
| -0.25 .. -0.20 | 5.03e-02 | 0.592 [0.561, 0.671] | 0 | 0.929 [0.370, 1.350] | 0 |
| -0.20 .. -0.15 | 6.43e-02 | 0.670 [0.606, 0.732] | 0 | 0.912 [0.450, 1.441] | 0 |
| -0.15 .. -0.10 | 7.57e-02 | 0.765 [0.684, 0.896] | 0 | 0.870 [0.620, 1.395] | 0 |
| -0.10 .. -0.05 | 8.84e-02 | 0.885 [0.761, 0.980] | 0 | 0.896 [0.730, 1.336] | 0 |
| -0.05 .. 0.00 | 1.01e-01 | 0.994 [0.908, 1.028] | 0 | 1.032 [0.917, 1.285] | 0 |
| 0.00 .. 0.05 | 1.39e-01 | 1.129 [1.059, 1.167] | 1.178 [1.048, 1.262] | 1.136 [1.070, 1.286] | 0.981 [0.857, 1.375] |
| 0.05 .. 0.10 | 1.92e-01 | 1.077 [1.017, 1.189] | 0.858 [0.755, 0.991] | 1.101 [1.035, 1.281] | 0.826 [0.716, 0.916] |
| 0.10 .. 0.15 | 2.60e-01 | 1.065 [1.020, 1.099] | 0.768 [0.712, 0.832] | 1.078 [1.032, 1.199] | 0.788 [0.593, 0.870] |
| 0.15 .. 0.20 | 3.39e-01 | 1.060 [1.037, 1.075] | 0.877 [0.831, 0.895] | 1.065 [1.034, 1.155] | 0.871 [0.755, 0.932] |
| 0.20 .. 0.25 | 3.99e-01 | 1.059 [1.037, 1.099] | 0.926 [0.889, 0.964] | 1.058 [1.031, 1.164] | 0.919 [0.857, 0.959] |
| 0.25 .. 0.30 | 4.62e-01 | 1.045 [1.031, 1.054] | 0.943 [0.922, 0.949] | 1.045 [1.025, 1.116] | 0.938 [0.888, 0.971] |
| 0.30 .. 0.35 | 5.24e-01 | 1.040 [1.026, 1.050] | 0.961 [0.942, 0.968] | 1.039 [1.021, 1.098] | 0.956 [0.919, 0.982] |
| 0.35 .. 0.40 | 5.98e-01 | 1.032 [1.019, 1.043] | 0.973 [0.956, 0.976] | 1.032 [1.017, 1.076] | 0.970 [0.942, 0.989] |
| 0.40 .. 0.45 | 6.94e-01 | 1.026 [1.016, 1.038] | 0.980 [0.972, 0.987] | 1.025 [1.013, 1.061] | 0.981 [0.960, 0.993] |

**Near the terminator (|mu0| < 0.15, n = 24 per side).** The median P_new/P_real is 1.012 for lon > 0
and 1.046 for lon < 0. The median P_old/P_real is 0.356 and 0.373. The morning and evening sides agree
to 3 %, so there is no strong east-west asymmetry at this time.

### Planet budget

Beam power absorbed by the whole planet, over all 6144 columns, in units of (1-A) F\* a_p^2:

- NEW: 5.385
- OLD: 5.106
- NEW / OLD = 1.055

**Where the rules differ.** Measured by column power, OLD and NEW differ by more than 1 % everywhere
with mu0 < 0.75: the curvature correction reaches that far. That region holds 58 % of the NEW total and
56 % of the OLD total. The terminator band holds much less:

| band | share of NEW total | share of OLD total |
|---|---|---|
| mu0 < 0.45 | 24.8 % | 21.6 % |
| mu0 < 0.20 | 7.7 % | 4.9 % |
| mu0 < 0.10 | 4.0 % | 1.9 % |
| mu0 < 0 | 2.0 % | 0 |

**Real rays in the sampled band.** The real-ray power in -0.55 <= mu0 < 0.45 is estimated by scaling
each bin's NEW power (all columns) with that bin's median P_real/P_new. In that band:

| | power | ratio to real |
|---|---|---|
| NEW | 1.337 | 1.017 |
| OLD | 1.104 | 0.840 |
| real (estimate) | 1.314 | 1 |

Columns with mu0 >= 0.45 are not ray-traced; they are assumed to have real = NEW. On that assumption,
planet total NEW / real = 1.004 and OLD / real = 0.952.

## Reading

1. **mu0 > 0.** Both rules take the target column's profile. NEW's residual is a systematic +3 to +6 %
   for mu0 = 0.1 to 0.45 and +13 % in the last dayside bin. That residual is the horizontal-inhomogeneity
   error alone, because gate 1 removes geometry. OLD adds its secant and clamp error on top:
   - At mu0 < 0.1, facsw = 10 caps the slant depth, and the /facsw factor sets the incident flux to
     0.1 F\* instead of mu0 F\*.
   - For 0.05 < mu0 < 0.2 OLD is 12 to 23 % low. For 0 to 0.05 it is 18 % high.

   On the dayside the geometry error is the larger of the two, and the new rule removes it.
2. **mu0 < 0.**
   - OLD is dark, so it is 100 % low.
   - NEW is 1 % low at the terminator and 30 to 50 % low for -0.45 < mu0 < -0.15: real grazing rays
     deliver about twice the power NEW computes. The per-cell scatter is wide (p10 to p90 is 0.2 to
     2.8), so the vertical placement of the heating differs too, not only its total.
   - Here the horizontal-inhomogeneity error dominates NEW's error. It is still smaller than OLD's,
     which is total.
   - The cause (which neighbouring columns make the chords more transparent) was not diagnosed.
   - The twilight ring holds 2 % of the planet's beam power, so the budget effect is small. Locally,
     NEW still under-heats the twilight upper atmosphere by up to a factor of 2.
3. **Budget.** NEW is within 0.4 % of the real rays for the planet total and within 1.7 % in the
   terminator band. OLD is 4.8 % and 16 % low respectively.

## Caveats

* One dump (rot 283.28), one MHD run.
* 8 columns per bin.
* The rays sample opacity from the nearest cell.
* Sampling stops once every chain has tau > 60 for both the real ray and the control.
* The EOS table is not the run's own: it lacks the condensation flag.
* The OLD rule's ghost column is a proxy, worth about 2e-4 of the flux (README_beam.md).
* icut is a transcription of the Rosseland blend, not read from the run.
