# he4_presn on the cubed sphere: the MLT SUB-GRID FLUX, arms D and E
(2026-09-17, viper, branch `he4-presn-global`, binary `../../build_gpu_rg/src/athena`)

`../arms/README.md` closed on this diagnosis: from the MLT-built ic the upper FeCZ drains
from `t = 0` (`L_rad,out/L` 0.94 -> 0.66 by 0.06 turnover, rho at 0.97 R -44 % in a third
of a turnover) because the two-stream sees a flux divergence where the 1-D structure had
15-40 % of `L` carried by *convection*, and resolved convection cannot spin up in that
time.  A radiative-equilibrium restart does not exist for this star (super-Eddington).

This directory tests the remaining option: **carry the convective flux from `t = 0` with
the pgen's own MLT sub-grid flux, then hand it over to the resolved flow.**  The inner
boundary is arm B's (`inner_bc = wall`, `rt_bottom_flux = true`,
`rad_flux_inner = 1.305278e15`), which `../arms/` settled on.

---

## What the MLT flux code actually does

All line numbers are `src/pgen/red_giant.cpp` at this commit.  The whole thing lives in
the pgen's user source term, behind `if (mlt_alpha_ > 0.0)` (**:3074**), and runs once per
RK stage with that stage's `bdt`.

**It is a radial face flux added to the ENERGY equation only.**  `rg_mlt_div`
(**:3531**) does

    u0(m,IEN,k,j,i) += bdt*(A_{i} F_{i} - A_{i+1} F_{i+1})/V     (curvilinear)

from `fconv_`, the per-face convective flux built by `rg_mlt_flux` (**:3415**).  Nothing
is added to the momentum or mass equations, and the EOS radiation taper / `rt_rad_force`
are untouched by it -- it is a pure heat-flux divergence.

**Two closures, selected by `problem/mlt_mean` (default TRUE).**

* *Per-column* (`mlt_mean = false`, **:3405-3471**): each column's own face gradient,
  `F = rho cp T sqrt(g delta) l^2 x^{3/2} / (4 sqrt2 Hp^{3/2})` with `l = alpha Hp`,
  `x = grad - grad_ad`, capped at `L/(4 pi r^2)`, then multiplied by the **tau taper**
  `w = clamp(log10(tau/tau_lo)/log10(tau_hi/tau_lo), 0, 1)` (`mlt_tau_lo/hi`), and by
  `mlt_rmax`.  A chi-limiter and a `0.1 e/dt` cap bound the explicit diffusion.
* *Shell mean* (`mlt_mean = true`, **:3108-3400**): the closure is evaluated **once per
  radial shell** on the horizontally averaged state (MPI-reduced, **:3172**), for the
  measured reasons at **:594-609** (per-column `x ~ 1e-5` deep down is round-off, and the
  3/2 power turns its sign noise into a face-to-face checkerboard of order `L`).

**The shell-mean closure is a DEFICIT closure, and that is the hand-over.**  At each face
(**:3269-3390**):

1. `F_req = L/(4 pi r^2)`.
2. Schwarzschild in the non-differentiating form `grad_rad = 3 kappa p F_req /
   (16 sigma g T^4) > grad_ad` decides whether the shell is convective (**:3341**).
3. The **deficit** `D = max(0, F_req - F_cond - F_res - F_2s)` (**:3347**), where
   - `F_cond` is the **conduction operator's own discrete face flux**, replicated term for
     term from the same `w0` (**:3183-3210**) -- not an estimate of it;
   - `F_2s` is the **two-stream's own net face flux** (`rt_face_flux()`, **:3155-3163**);
   - `F_res` is the **RESOLVED convective flux**
     `<rho v_r h> - <h><rho v_r> + <rho v_r v^2/2>`, formed per cell and decorrelated
     (**:3292-3305**).
4. The applied target is `min(F_MLT(x), D)` where `x > mlt_x_thr` (the amplitude is
   trusted), and `D` alone where `x` is jitter (**:3353-3366**).
5. `mlt_relax_time` (default 1e4 s) relaxes the stored 1-D profile `fmlt1d(i)` toward
   that target by `dt/relax_time` per call; it is **seeded** on the first call, so `t = 0`
   already carries the full closure (**:3374-3381**).

So: **it is recomputed from the live state every step, not frozen from the ic**, and it
**dies on its own** as the resolved flow takes the load, because `F_res` is subtracted.
The tau taper and `mlt_rmax` are **not used** in this path -- "the deficit IS the
hand-over" (**:3241-3243**); `mlt_tau_lo/hi` are dead parameters under `mlt_mean = true`,
and `taumlt_` is not even filled (**:3094**, `if (!mltmean)`).  Per column only the safety
cap `0.1 min(e_l, e_r) dx / bdt` is re-applied (**:3415-3425**).

**The ramp went UP only.**  `mlt_ramp_time` multiplies the applied flux by
`min(1, t/mlt_ramp_time)` (**:565-570**) -- a build-up for a switch-on transient, the
opposite of a hand-over.

**Switch added (this commit).**  `problem/mlt_ramp_down_time` and `problem/mlt_hold_time`:
the applied flux is multiplied by 1 for `t <= mlt_hold_time` and falls linearly to 0 at
`mlt_hold_time + mlt_ramp_down_time`.  Default 0 = off, so every existing run is bitwise
unchanged; `red_giant.cpp` only.  The shell-mean closure already hands over *physically*,
but that hand-over is a fixed point, not a controlled one: if the resolved flow stalls the
subgrid flux simply carries `L` for ever and the run never tests whether convection
started.  The ramp-down forces the question.  With it on, rank 0 prints
`### red_giant: MLT subgrid flux fraction = <ramp> at t = <t> s` every 100 cycles.

**Interaction with the seed.**  `problem/vpert_mlt` defaults TRUE and is gated on
`mlt_alpha_ic`, not `mlt_alpha` (**:2569**), and the input already carries
`mlt_alpha_ic = 1.5`.  So arms A/B were already seeded with `vpert * v_c(r)`, and turning
`mlt_alpha` on **does not change the seed**.

---

## The runs

`chain.sh` (`--export=ALL,ARM=<D|E>,WALL=hh:mm:ss,TLIM=<s>`), `-p apudev`, 2 MI300A, the same
smoke grid as `../arms/` (nx1 = 96, nx2 = nx3 = 32, `meshblock 96x16x16`, 24 blocks).
Tables from `../arms/gate.py <armdir>`.

| arm | `mlt_alpha` | `vpert` | hand-over | tlim |
| --- | --- | --- | --- | --- |
| D | 1.5 | 1e-3 | none (held) | 1 turnover = 4705 s |
| E | 1.5 | 1e-2 | hold 1 turnover, ramp to 0 over turnovers 1-4 | 5 turnovers = 2.3525e4 s |

## ARM D -- MLT flux held, `vpert = 1e-3`, one turnover -- **FAIL**

Job 11759591, `apudev` (2 MI300A on vipa1327; the `apu` partition went into maintenance
at 10:00).  `tableD.txt`, generated by `../arms/gate.py D`:

```
  t/turn      t[s]     dt[s]    lnKE1   lnKEh   dmass    eos_fail floors  fofc   L/Lstar   d_rho_0 d_T_0  d_rho_1 d_T_1  d_rho_2 d_T_2   vr_rms/vMLT
   0.00         0      12.9   75.50   74.64  0.0e+00        0       0      0     0.955   0.000   0.000   0.000   0.000   0.000   0.000     0.000
   0.04     192.8      14.4   87.23   74.64  1.4e-07        0       0      0     0.742  -0.124   0.022  -0.009  -0.023   0.095  -0.045     0.003
   0.08     382.1      12.9   88.50   76.27  2.4e-07        0       0      0     0.673  -0.001   0.033  -0.092  -0.024   0.182  -0.056     0.159
   0.12     564.8      11.6   88.72   79.72  2.2e-07        0       0      0     0.688   0.486   0.039  -0.119  -0.024   0.174  -0.061     0.270
   0.16     755.1      10.7   88.69   80.22  2.8e-08        0       5      0     0.700   0.256   0.043  -0.151  -0.023   0.130  -0.068     0.189
   0.20     945.7      10.7   88.52   81.62 -1.6e-07        0       8      0     0.714   0.238   0.046  -0.160  -0.022   0.071  -0.080     0.054
   0.24      1130      10.9   88.91   82.95 -3.1e-07        0      12      0     0.721   0.234   0.048  -0.157  -0.019   0.008  -0.094     0.081
   0.28      1326      10.8   89.06   83.47 -4.2e-07        0      71      0     0.730   0.451   0.047  -0.155  -0.015  -0.066  -0.112     0.154
   0.32      1504      8.28   89.39   84.13 -5.7e-07        0     299      0     0.761   0.375   0.045  -0.163  -0.012  -0.138  -0.127     0.188
   0.36      1694      6.29   89.57   85.50 -7.3e-07        0     720      0     0.779   0.230   0.045  -0.165  -0.008  -0.214  -0.147     0.220
event-log totals: eos_fail=0 dfloor=0 efloor=720 tfloor=0 vceil=0 fofc=0 max c2p it=0
```

`dt` collapses at **cycle 164, t = 1723.65 s = 0.366 turnover** -- arm B died at
**t = 1723.03 s**.  The first collapsing cell is `(m,k,j,i) = (7,18,9,22) gid = 19`,
`r = 1.593e11 = 0.672 R`, `rho = 1e-13` (ON `dfloor`), `p = 0.1` (ON `pfloor`),
`v1 = -4.1e22`; arm B's was `i = 24`, `r = 1.638e11 = 0.690 R`.  **Same shell, same
second.**

### The gate

| gate | required | measured | verdict |
| --- | --- | --- | --- |
| `(L_rad + F_MLT)_out / L` within 5 % of 1 | 0.95-1.05 | 0.955 at t = 0, **0.673** by 0.08 turnover, 0.78 at the death | **FAIL** |
| `efloor` = 0 through one turnover | 0 | 0 to t = 570 s, then 5, 12, 71, 299, **720** by 0.36 turnover | **FAIL** |
| rho / T at the FeCZ top drifting < 10 % | < 0.10 | rho **-21.4 %**, T **-14.7 %** at 0.36 turnover | **FAIL** |
| reaches one turnover | 4705 s | 1723 s | **FAIL** |

The MLT flux is **not neutral** -- it measurably helps, it just does not help enough:

| at 0.36 turnover | B (no MLT) | D (MLT held) |
| --- | --- | --- |
| `L_rad,out/L` | 0.721 | **0.779** |
| rho at 0.97 R | -43.6 % | **-21.4 %** |
| T at 0.97 R | -18.2 % | **-14.7 %** |
| `efloor` | 1230 | **720** |
| `lnKEh` | 83.59 | **85.50** (6.7x the horizontal KE) |
| death | t = 1723.03 s | t = 1723.65 s |

`eos_fail`, `dfloor`, `tfloor`, `vceil`, `fofc` and `max c2p it` are all **exactly zero**
for the whole run, as in arm B.  Mass drift is -7.3e-7.

### Which term -- the `mlt_dump` face budget at t = 0

This is what `problem/mlt_dump` was for.  `D/mltfaces_D.txt`, every sixth face, all fluxes
as fractions of `F_req = L/(4 pi r^2)`.  `F_raddiff` is the plain Rosseland diffusion flux
`16 sigma T^3/(3 kappa rho) |dT/dr|` on the **same shell-mean face state**; `F_2s` is the
grey two-stream's **own** net face flux; `F_MLT` is the mixing-length amplitude and
`F_used` what the closure applied.

| i | r/R | x = grad-grad_ad | F_raddiff/F_req | F_2s/F_req | F_MLT/F_req | F_used/F_req | (F_2s+F_used)/F_req |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 4 | 0.506 | -1.08e-02 | 1.071 | 0.735 | 0.000 | 0.000 | 0.735 |
| 10 | 0.553 | -6.09e-03 | 1.017 | 0.481 | 0.000 | 0.000 | 0.481 |
| 16 | 0.609 | -3.33e-03 | 1.050 | 0.182 | 0.000 | 0.000 | 0.182 |
| 22 | 0.667 | +9.12e-04 | 0.986 | 0.107 | 0.050 | 0.050 | 0.157 |
| 28 | 0.722 | +1.95e-03 | 0.910 | 0.160 | 0.120 | 0.120 | 0.280 |
| 34 | 0.772 | +3.18e-03 | 0.867 | 0.241 | 0.171 | 0.171 | 0.413 |
| 40 | 0.815 | +4.80e-03 | 0.851 | 0.336 | 0.203 | 0.203 | 0.539 |
| 46 | 0.851 | +7.18e-03 | 0.861 | 0.440 | 0.216 | 0.216 | 0.656 |
| 52 | 0.880 | +1.13e-02 | 0.888 | 0.551 | 0.223 | 0.223 | 0.775 |
| 58 | 0.902 | +1.94e-02 | 0.929 | 0.681 | 0.223 | 0.223 | 0.905 |
| 64 | 0.921 | +3.66e-02 | 0.961 | 0.838 | 0.180 | 0.162 | 1.000 |
| 70 | 0.937 | +3.71e-02 | 0.976 | 0.929 | 0.038 | 0.038 | 0.967 |
| 76 | 0.952 | +1.82e-02 | 0.986 | 0.957 | 0.003 | 0.003 | 0.960 |
| 82 | 0.968 | +7.24e-04 | 0.990 | 0.969 | 0.000 | 0.000 | 0.969 |
| 88 | 0.984 | -3.38e-03 | 0.972 | 0.983 | 0.000 | 0.017 | 1.000 |
| 94 | 1.001 | -1.44e-02 | 0.797 | 1.185 | 0.000 | 0.000 | 1.185 |

(`F_cond/F_req = 0` and the blend weight `w = 0` at **every** face: `rad_tau_lo/hi =
1e5/1e6` hands the whole column to the two-stream, so the radial diffusion operator is
inert by construction.)

Read the last column.  **At t = 0, between 0.55 R and 0.92 R, the run carries only
16-90 % of the luminosity the star has to carry, with a minimum of 0.157 at r/R = 0.667**
-- which is where both arm B and arm D put their first collapsing cell.  The flux rises
monotonically outward through that range, i.e. the shells there are being drained, which
is exactly the measured behaviour (upper FeCZ empties, `L_rad,out/L` falls to 0.67, the
cell at 0.67 R ends on `dfloor` and `pfloor`).

**The term is the two-stream's own net radial flux in the optically thick interior, and
it is NOT the initial structure and NOT the missing convection.**

* `F_raddiff/F_req` is **0.85-1.07 at every interior face**.  Radiative diffusion, with
  this run's own opacity table, on this run's own shell-mean state, carries essentially
  the whole luminosity everywhere below the photosphere.  The ic is self-consistent.
* `F_2s/F_req` on the identical state is **0.107-0.74** over 0.55-0.90 R, a factor
  **2-9 short**, and only reaches 0.96 above 0.95 R.  In the diffusion limit the two must
  agree: the per-cell optical depth at the base is `kappa rho dr = 0.62 x 2.36e-8 x 6e8
  ~ 9`, and the whole domain is only `tau = 349` deep, so this is a regime where the grey
  two-stream should reproduce `(4 pi/3) dB/dtau` exactly.
* Consequently the sub-grid model was never the missing carrier.  Where it is allowed to
  act (`grad_rad > grad_ad`, i.e. `r/R > 0.63`) it supplies 0.05-0.22 `F_req` -- the
  MLT **amplitude** binds, not the deficit, because the shell-mean superadiabaticity of
  this structure is only `x = 1e-3 .. 4e-2`.  And below 0.63 R the shells are
  Schwarzschild-**stable** (`grad_rad < grad_ad`, `x < 0`), so the closure correctly
  refuses to supply anything at all -- yet that is exactly where the budget is worst
  (0.18 `F_req` at 0.609 R).  **No convective closure, sub-grid or resolved, can fix a
  hole in a stable layer.**

## ARM E -- NOT RUN

Its precondition was arm D passing.  Per the brief, the work stops at the diagnosis.
The `mlt_ramp_down_time` / `mlt_hold_time` switch is in and tested only in the sense that
arm D ran with it off (default, bitwise unchanged).

## Verdict

**The MLT sub-grid flux does not close the budget, and it cannot: the hole is in a
Schwarzschild-stable layer where the grey two-stream, not convection, is supposed to carry
the flux, and it carries 2-9x too little.**  Arm D dies within 0.6 s of arm B.

## Recommended next step

**Stay on the smoke grid and test the radial transport operator, not the boundary, not the
seed, and not the convection model.**  In order, cheapest first (each is one 15-minute
`apudev` job, and the first two need no code change):

1. **Give the deep interior back to radiative diffusion.**  `rad_tau_lo/hi = 1e5/1e6` was
   copied from the plane-parallel He box to make the two-stream own all 96 cells; put the
   handover back at a sane optical depth (e.g. `rad_tau_lo/hi = 20/300`, with
   `rad_blend_radial`) so the conduction operator -- whose flux is `F_raddiff`, i.e. the
   one that IS right -- carries the interior, and the two-stream owns only the thin
   layers.  Gate on the same `mlt_dump` budget at t = 0: `(F_cond + F_2s)/F_req` must be
   within a few per cent of 1 at every face **before** any time-stepping.  This one table
   is the whole test and costs one cycle.
2. **Isolate the discrepancy inside the two-stream.**  With the budget dump as the meter,
   vary `rt_impl_mixed` (2 -> 0, single-precision factors off) and `ck_nquad` (2 -> 4),
   and compare `F_2s` against `F_raddiff` face by face.  If neither moves it, the
   disagreement is in the grey column solver's diffusion limit and belongs in a 1-D unit
   test (`tests_1d`) against the analytic `tau = 349` grey slab, not in this star.
3. Only when `(F_carried)/F_req = 1 +- few %` at t = 0 does any of the previous work
   become testable again: re-run arm B, then arm D, then arm E (the ramp-down) on this
   same grid.  Nothing about the production grid should be attempted before that.

## Files

* `chain.sh` -- the submit script (`--export=ALL,ARM=<D|E>,WALL=hh:mm:ss,TLIM=<s>`).
* `D/` -- `he4.hydro.hst`, `he4.log` (event log), `rt_profile.bin`, `column_he4_D.txt`,
  **`mltfaces_D.txt`** (the t = 0 face budget above).  `.bin`/`.rst`/`.cbin`/
  `rt_surface.bin` deleted after measuring.
* `tableD.txt` -- the gate table as generated.
* `D.11759591.log` -- the job log, truncated (startup block, cycle diagnostics, the
  dt-collapse block); every number above was read off the full log first.

