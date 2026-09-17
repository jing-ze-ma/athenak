# The SPHERICAL FORM of the grey two-stream: gates G1-G5

2026-09-17, viper, branch `he4-presn-global`, commit `8bca3dfa` against HEAD `75b95308`.
Binaries: `build_cpu_box` (serial CPU, `PROBLEM=box_convection`, G1) and `build_gpu_rg`
(MI300A, `PROBLEM=red_giant`, G2-G5).  References were built from `75b95308` first and
kept; both have been deleted here after measuring, and so have every `.bin`, `.rst` and
`.cbin`.  Every number below was read off the full outputs before the pruning.

## What changed, in one paragraph

`1159a8f3` transported `J = A I` against a source `A B`.  That puts the area on BOTH
two-stream moments; the exact spherically symmetric pair with `f = mu^2 = 1/3` puts it on
the antisymmetric one alone,

    mu d(A D)/dr = -kappa rho A (S - B),      mu dS/dr = -kappa rho D,

so the old scheme's diffusion flux was `F_2s/F_exact = 1 - H_T/(2r)`.  Now each cell
carries a constant frame area `A_i = V_i/dx_i`, the layer solve inside it is
plane-parallel per unit area with source `B`, and at every face `S` and `A D` are passed
continuously -- one correction `c = beta (d_above - u_below)` added to both sides, with
`beta = (A_above - A_below)/(A_above + A_below)`.  Full derivation: the SPHERICAL FORM
note on `RTCol3::Bt` in `src/utils/two_stream_column_implicit.hpp`.

## G1 -- the plane-parallel box regression: BITWISE

`bench/hestar_fecz/box_w8/he_box_w8.athinput`, production mode-3 configuration
(`rt_implicit_column = 3`, `rt_impl_solver = pcr`, `rt_impl_mixed = 2`,
`rt_impl_warm = 1`, `rt_col3_skip_sweep = true`), shrunk to `nx2 = nx3 = 16`,
`meshblock 134x8x8`, **50 cycles**, serial CPU.  `g1.sh` is the driver.

    cmp g1_ref/feczrt.hydro.hst  g1_new/...   IDENTICAL
    cmp g1_ref/feczrt.user.hst   g1_new/...   IDENTICAL
    cmp g1_ref/column_used.txt   g1_new/...   IDENTICAL
    cmp g1_ref/rt_surface.bin    g1_new/...   IDENTICAL
    cmp g1_ref/rt_profile.bin    g1_new/...   IDENTICAL
    diff -r bin/ cbin_hydro_w_2/ rst/         IDENTICAL

and the same six lines again for the **mode-0 variant** (`g1m0_ref` / `g1m0_new`,
`rt_implicit_column = 0`, `rt_col3_skip_sweep = false`).  Byte-identical, not
"agrees to round-off": every new area factor is a multiplication by a `1.0` READ FROM THE
AREA VIEW (`Aun`/`Acn`/`AUN`/`ACN`), never a folded literal and never a division by an
exact 1.0, and `beta` is exactly 0, so no expression is reshaped and hipcc's
`-ffp-contract=fast` sees the same DAG.

## G2 -- the transparent gate: PASS

`tests_m3/thin_ex.athinput` (= `inputs/hydro/red_giant_1d_dilution.athinput` +
`rt_col3_ex_iter`): cubed sphere, `nx2 = nx3 = 4`, 128 uniform radial cells,
`r_out/r_in = 2.5` so the area ratio is 6.25, `kappa_const = 1e-5`, `vpert = 0`.
`L = 4 pi r^2 F` from the solver's own face flux, read with `tests_r2/rg1d/lprof.py`.
Gate: `max/min < 1.005`.

| arm | `L` max/min |
| --- | --- |
| mode 0 | **1.00031** |
| mode 3 `thomas`, `rt_col3_skip_sweep` | **1.00018** |
| mode 3 `pcr`, `rt_col3_skip_sweep` | **1.00018** (identical to `thomas` in every printed digit) |

## G3 -- the thick gate: PASS, and it is the one the old scheme failed

`inputs/tests/two_stream_sph_thick.athinput`, all eight cases of
`tests_r2/thick/cases.txt`, one cycle each.  `F_2s/F_req` over the interior faces with
`dtau_cell > 0.3`, the wall face and the top ghost mirror excluded exactly as
`tests_r2/thick/README.md` section 3 does.  `an3b.py` builds `tab_g3.txt`:

| case | mode 3: max\|1-x\| / median | mode 0: max\|1-x\| / median | HEAD `J = A I` |
| --- | --- | --- | --- |
| n0_t100_r2.0 | 0.0011 / 0.99995 | 0.0007 / 0.99995 | 0.972 / 0.461 |
| n1_t100_r2.0 | 0.0047 / 0.99994 | 0.0032 / 0.99995 | 0.738 / 0.531 |
| n3_t100_r2.0 | 0.0020 / 0.99990 | 0.0046 / 0.99993 | 0.467 / 0.627 |
| n7_t100_r2.0 | 0.0004 / 0.99974 | 0.0004 / 0.99975 | 0.256 / 0.748 |
| n3_t30_r2.0  | 0.0001 / 0.99989 | 0.0001 / 0.99999 | 0.473 / 0.538 |
| n3_t300_r2.0 | 0.0033 / 0.99990 | 0.0004 / 0.99991 | 0.465 / 0.642 |
| n3_t100_r1.1 | 0.0025 / 0.99995 | 0.0014 / 0.99995 | 0.158 / 0.907 |
| n3_t100_r1.5 | 0.0120 / 0.99994 | 0.0083 / 0.99994 | 0.399 / 0.721 |

The gate was 1 %.  The median is 1e-4 from 1 in every case and the worst single face is
1.2 %; the old scheme was 16-97 % off, exactly the `1 - H_T/(2r)` the directory measured.
**Mode 0's lagged approximation costs nothing measurable here** -- its worst case (0.83 %)
is no worse than mode 3's (1.2 %), because `beta = dr/2r ~ 4e-3` and the lag is `O(beta^2)`.

**A trap this exposed.** `F_2s` (column 17 of `mlt_dump`) is `Fb`, which the EXPLICIT
sweep writes unless `rt_col3_skip_sweep` tells the column solve to own it.  The first
measurement of `tests_r2/thick` -- and its README's "MEASURED" tables -- therefore report
**mode 0's** flux whatever `rt_implicit_column` says.  The input now carries
`rt_col3_skip_sweep`, `rt_col3_ex_iter` and `rt_src_direct` with a note, and the mode-3
arm above passes `problem/rt_col3_skip_sweep=true`.

## G4 -- the He star `t = 0` face budget: PASS

`inputs/hydro/he4_presn_cs.athinput` on the smoke grid AS COMMITTED (whole-column blend,
`rad_tau_lo/hi = 1e5/1e6`, so `w = 0` on every face and the two-stream is the only
carrier), `inner_bc = wall`, `rt_bottom_flux = true`,
`rad_flux_inner = 1.305278e15 = L/(4 pi r_in^2)`, `vpert = 0`, `mlt_alpha = 1.5`
(which is what makes `mlt_dump` fire), `time/nlim = 1`.  `an4.py` builds `tab_g4.txt`.
The two INDEPENDENT carriers of the same flux on the SAME state, over 0.55-0.97 R
(73 faces):

| | `F_2s/F_req` | `F_raddiff/F_req` | `F_2s/F_raddiff` |
| --- | --- | --- | --- |
| HEAD `75b95308` | 0.102 .. 0.969 (med 0.457) | 0.851 .. 1.059 (med 0.949) | **0.099 .. 0.979** (med 0.511) |
| this commit | 0.843 .. 1.003 (med 0.940) | 0.850 .. 1.008 (med 0.942) | **0.983 .. 1.016** (med 0.998) |

Gate was 5 %.  The 2-9x shortfall of `tests_3d/handover/README.md` section 1 is gone: the
grey two-stream and Rosseland diffusion now agree to 1.7 % at every face of the envelope.

## G5 -- the 3-D smoke arms

`chain5.sh`, the `../tests_3d/arms` grid and cadences: `nx1 = 96` stretched from 0.50 R,
`nx2 = nx3 = 32` per panel, `meshblock 96x16x16` = 24 blocks, 2 MI300A, turnover
= 4705 s.  Both arms are the old **arm B** boundary (`inner_bc = wall`,
`rt_bottom_flux = true`, `rad_flux_inner = 1.305278e15`) with `mlt_alpha = 0`.
`../tests_3d/arms/gate.py` builds `g5/tableBp.txt`, `g5/tableEp.txt`.

### B' -- `vpert = 1e-3`, one turnover: **SURVIVES** (the old arm B died at 0.366)

| t/turn | dt[s] | eos_fail | floors | fofc | `L_out/L` | d_rho base | d_rho peak | d_rho FeCZ-top | `v_r rms / v_MLT` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.00 | 12.9 | 0 | 0 | 0 | 0.996 | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.20 | 12.8 | 0 | 2.4e5 | 0 | 0.978 | +0.001 | -0.022 | -0.030 | 0.078 |
| **0.36** | 12.2 | 0 | 2.5e6 | 0 | **0.989** | **0.000** | **-0.072** | **-0.091** | 0.136 |
| 0.50 | 11.2 | 0 | 3.3e6 | 0 | 0.988 | -0.001 | -0.082 | -0.158 | 0.099 |
| 0.70 | 9.69 | 0 | 3.6e6 | 0 | 0.974 | -0.002 | -0.149 | -0.277 | 0.150 |
| 1.00 | 7.1 | 0 | 3.6e6 | 0 | **0.917** | -0.003 | -0.292 | **-0.544** | 0.123 |

`eos_fail = 0`, `fofc = 0`, `tfloor = 0`, `vceil = 0`; `dfloor = 2.29e6`,
`efloor = 1.30e6` cell-cycles.  Run reached `tlim` and exited 0.
**The old arm B at the same time (0.36 turnover, where it died):** `L_out/L = 0.721`,
`d_rho` base **+0.280**, peak +0.150, FeCZ-top **-0.436**.  So at equal time the new run
is better by 4-30x on every structural number, and it then keeps going.

Against the gate as written: alive **yes**; `eos_fail = 0` **yes**; `L_out/L` within 10 %
**yes** (0.917-0.996); `efloor = 0` **no** (1.3e6 cell-cycles, against 1230 for the old
arm at its death -- but the new run is 3x longer and its floors sit in the draining
top atmosphere, not at a collapsing cell); FeCZ-top `rho` drift < 20 % **no** at one
turnover (-54 %), **yes** at the 0.36 turnover the old arm reached (-9.1 %).

### E' -- `vpert = 1e-2`, three turnovers: dies at **0.782** turnover

Mean structure tracks B' face for face (same `L_out/L`, same `d_rho`, same floors) -- only
`lnKE_h` differs, by exactly the 100x the seed puts in.  First collapse:

    ### dt COLLAPSE cycle=328 time=3659.33 dtold=4.40747 dt=1.04179 | hydro=1.04179
        hydro dt is set by cell (m,k,j,i) = (8,18,18,30) gid = 20
        r=1.76452e+11 (r/R = 0.744, the opacity peak) rho=3.52449e-09 T=5.24607e+11
        p=459.013 cs=463942 v=(5.76332e+08,-65690.6,-37142.5)

`p = 459` at `rho = 3.5e-9` and `T = 5.2e11 K` is not a thermodynamic state: ideal `p`
there would be 1e8.  It is the EOS inversion on a floored internal energy -- the
`red-giant-seam-floor-eos-garbage` signature -- so the trigger is the **floor**, in the
same upper FeCZ that B' is slowly draining, and not the transport.  The face budget at
that moment is still healthy (`L_rad,out/L = 0.967`, `L_rad,cut/L = 1.006`); it is
`-3.1e5` sixteen seconds later.

**Conclusion of G5.** The transport defect is fixed and it was worth 3x in survival time
and 30x in the base-density error, but it is not the only thing wrong with this star: the
upper FeCZ drains, hits `efloor`, and a strong enough seed turns that into a detonation.
The next thing to chase is the floor in the draining layer, not the radiation operator.

## Files

* `g1.sh`, `g1_ref/ g1_new/ g1m0_ref/ g1m0_new/` -- G1, the `.hst` and event logs kept.
* `g2/{m0,m3s,m3p}/run.log` -- G2 (the `L` profile is in the log).
* `g3_n3/ g3_n0/ g3_r0/` -- G3, `mltfaces.txt` only; `an3b.py` -> `tab_g3.txt`.
* `g4_sw/ g4_ref/` -- G4, `mltfaces.txt` + `he4.log`; `an4.py` -> `tab_g4.txt`.
* `chain5.sh`, `g5/Bp/ g5/Ep/` -- G5, `he4.hydro.hst`, the event log, `rt_profile.bin`,
  the truncated job logs, and `tableBp.txt` / `tableEp.txt`.
