# HANDOVER 2026-09-17 B (viper, evening): the global 4 Msun presupernova He star — code done,
# a spherical two-stream defect blocks the star

Read after `HANDOVER-2026-09-17.md` (the 02:00-03:00 one) in this same directory.  This document
covers 2026-09-17, 07:00-18:40 CEST, on viper only.
Memory snapshot: `docs/handover/claude-memory-2026-09-17-B/` (copy into the local memory dir;
`MEMORY.md` first).  Analysis scripts from earlier days: `docs/handover/scripts/`.
Runs live in `/viper/u2/jinma/ATHENAK/bench/` (viper-only).  Never write into `run/`.

**One sentence:** the whole global-sphere pipeline for the 4 Msun presupernova He star was built,
merged and gated today (36 commits on `he4-presn-global`), and then every 3-D smoke arm died inside
0.4 turnover on a single root cause that is now measured exactly — the grey two-stream, made
spherical by transporting `J = A I`, is **wrong in the diffusion limit by `1 - H_T/(2r)`**, an O(1)
factor wherever `H_T ~ r`.  The fix is a moment-equation rewrite with a variable Eddington factor;
two unit tests to gate it already exist.  **Nothing was pushed.**

## 1. The two productions (both healthy, both untouched today)

| | B star `bench/bstar_fecz/prod_w7` | He star `bench/hestar_fecz/box_w8` |
| --- | --- | --- |
| at 18:38 | `t = 3.7270e5` s = **37.3 / 40 turnovers** (1 turnover = 1e4 s) | `t = 9.7184e4` s = **25.8 / 40 turnovers** (1 turnover = 3.761e3 s) |
| running job | `11735625` `feczrt_pw7n4`, 4 nodes, 1:53 in, `11735626` queued behind it | `11740424` `hefeczrt_bw8n4`, 4 nodes, 1:53 in, `11740425` queued |
| dumps | `bin` 37 | `bin` 25 |
| state | developed convection since ~turnover 15; 23-turnover check: rms v_z/v_MLT 0.73 / 0.88 at tau 300 / 1000 | 13-turnover check: rms v_z/v_MLT 0.93 (tau 20), 0.90 (tau 50), closed network of 4-5 cells |
| disk | **56 GB**, 28 `rst` | 8.2 GB, 48 `rst` |

`box_w9` (the stronger He seed) was **cancelled and deleted** by the user this morning: `box_w8`'s
growth turned out to be accelerating, not slow (ln-KE ratios over turnovers 8-13: +0.64 +0.25 +0.02
+0.09 +0.24 +0.18), so w9 was redundant.  `box_w8` hung at 08:09 on a zero-byte `cbin` write during
the quota event below; it was cancelled at 08:45 and the chain restarted from `rst 00026`.

Page **v7 version 3** https://claude.ai/artifact/SLgYFj9JDBUsNqxj5UYgPJ carries B star dump 23 and He
w8 dump 13 (plus the earlier 3, 10, 15 and 7); the 3-D viewers were deliberately left on older data (user).  Source `athenak/.pages_w7/`.

**Standing rule from today (user):** judge a convection run against **v_MLT and F_conv/F**, never
against an earlier run.  `box_w4` was 19x v_MLT (over-driven by the old whole-column kick), so the
earlier "200x below w4" verdicts meant nothing.  See `compare-with-mlt-not-old-runs.md`.

## 2. What was built for the global He star

Target (Woosley 2019, 4.0 Msun presupernova): `M = 3.15 Msun`, `log L = 4.78`, `Teff = 49 kK`,
`R = 2.3717e11 cm`, `L = 2.3066e38`.  Sized in `bench/hestar_presn/GAP_ANALYSIS.md`: **no**
presupernova He star is a box problem — the FeCZ spans 0.42-0.74 R, `Hp/R = 0.3-0.7`, `Prad/Pgas`
27-160, `kappa F/(c g) = 1.1-2.4` (super-Eddington), so only a global spherical RHD model works.

Design decisions, all taken today:

- **Base = `src/pgen/red_giant.cpp`**, not `box_convection.cpp`: port the box physics (~1000-1200
  lines) onto the existing cubed-sphere machinery rather than teach the box spheres.  Envelope mass
  is 4e-5 Msun, so the point mass is exact.
- **Full sphere, cubed sphere**, `use_cubed_sphere = true`, one block per panel radially (whole
  column in one MeshBlock — the mode-3 column solve requires it).
- **2 cells per surface Hp horizontally** -> `nx2 = nx3 = 320` per panel (`Hp_surf = 3.086e9`,
  largest gnomonic cell `2r/N`); production `nx1 = 192` stretched.
- **Inner boundary at 0.5 R** (user, 08:30; `GAP_ANALYSIS` had 0.4 R): 0.4-0.5 R is 6.6 Hp of stable
  0.24-1 MK gas and would eat half the radial cells; 1 Hp / 2 H_rho of buffer remain below the FeCZ
  base.
- **Strang splitting** for the radiative operator (`rt_strang = true`), mode-3 exact column solve
  (`rt_implicit_column = 3`, `pcr`, `rt_impl_mixed 2`), `rt_rad_force`, FOFC, top + bottom tau-based
  sponges placed clear of the FeCZ (FeCZ `tau` 11-208 on the sphere), implicit cs ADI transverse
  operator.
- **Pipeline** (`bench/hestar_presn/`): `column_sph.py` -> `column_he4_presn_sph.txt` (spherical 1-D
  hydrostatic + MLT march) -> `make_ic_sph.py` -> `ic_he4_presn_sph.txt` (r, rho, eint with the
  taper applied, 0.35-1.026 R) -> the new `problem/ic_profile` reader in `red_giant.cpp`.  This
  route exists because red_giant's **built-in** 1-D march cannot build a Prad-dominated envelope
  (fatal `aT^4/3 > 0.9 p_top`, ignores taper/force, 1e6 rho error at `r_in`).
- Spherical structure: FeCZ **0.635-0.968 R** (tau 11-208), `v_MLT = 1.45e7` cm/s (Ma 0.42),
  `F_conv/F = 0.16`, **turnover 4705 s**, thermal time 5.7 d, a physical **density inversion** over
  0.72-0.95 R (contrast 3.9, present in both columns).  The plane-parallel column is only valid
  above 0.98 R.
- Diagnostics added: `rt_surface` (panel, x2v, x3v, F_top rows), `rt_profile` (8 shell means),
  `mlt_dump` (per-face flux budget), `column_dump`.

Gates that **passed**: `ic_profile` read-back 6.3e-9 in rho, 6.5e-7 in p; the run's own
`tau = 2/3 at r = 2.37122e11, T = 49199 K` against the model's `R = 2.3717e11, Teff = 48978 K`
(0.02 % / 0.45 %); red_giant bitwise unchanged by all 56 new switches; box production config bitwise
through every change.

Grid and cost for the first production run (GAP_ANALYSIS §E, recalibrated on `prod_w7`):
`6 x 320^2 x 192 = 1.18e8` cells, `dt ~ 3 s`, ~2000 cycles/turnover, `1.6e7` zone-cycles/s/GPU
measured with the full physics -> **~1 h wall, 8 GPU-h per turnover** on 8 MI300A with a 2x
cubed-sphere penalty.  Smoke grid: `nx1 = 96`, `nx2 = nx3 = 32`, blocks `96x16x16` = 24, 2 MI300A,
under a minute per turnover.  Since the seam fix (§3) **2x2 blocks per panel is no longer the limit —
4x4 is now clean**, so the production block decomposition can be chosen freely.

## 3. Bugs found and fixed today

| commit | bug | who is affected |
| --- | --- | --- |
| `1159a8f3` | the two-stream deposited `-(F_top-F_bot)/dx1` with **no `r^2` dilution** (factor 6.4 over `r_out/r_in = 2.54`) | every spherical two-stream run; box bitwise.  Fixed by area-weighted intensities `J = A I` + volume deposit + mode-3 rows.  1-D spherical `L` const to 0.03 % vs x6.25 before.  **This is also the commit that introduced the blocker of §4.** |
| `bbcb06a0` | `1159a8f3` carried **absolute** areas (`A ~ 1e22` next to unit entries) into the mode-3 5x5 blocks -> the single-precision pivot test failed on every block and mixed precision was silently disabled; on the sphere mode 3 blew up (`L_out/L = 8.6e5`) | fixed by carrying areas **relative to the column top**.  Radial arithmetic of the mode-3 residual is otherwise verified correct to round-off (`tests_m3/`). |
| `4bcdc855` | **`rt_rad_force` transverse term**: `Prad grad w` divided by the **coordinate** spacing `size.dx2/dx3` instead of the arc length.  On the cubed sphere x2 is the gnomonic angle in [-1,1] -> the force was **2e11x too large**; on **spherical polar** `dx2 = dtheta` -> too large by a factor **r** | fixed with `DX2/DX3` arc-length helpers (Cartesian bitwise).  **Every `red_giant` / `deep_hot_jupiter` run on spherical polar with `rt_rad_force = true` carried this.**  Check whether the dhj/rg productions used it before trusting their thin-layer horizontal dynamics.  **Not yet ported to `rt-integration` / `polar-average-perf`.** |
| `94c7165d` | the cubed-sphere **cell-centred seam halo was wrong with >1 MeshBlock per panel** (2x2: 11x, 4x4: 100x worse): (1) the monotonicity clamp was applied unconditionally to an *extrapolating* stencil at interior block ends, (2) the x2x3 edge ghosts (slots 40-47) across a seam were a plain copy, never resampled | after the fix 2x2/3x3/4x4 = 1x1 to printed digits, 1x1 bitwise, 4 cs tests pass, blast dt history unchanged.  **`src/bvals/bvals_fc.cpp` (face-centred B) has BOTH defects and is untouched -> every cs MHD run with >1 block per panel, including the dhj cs MHD productions, carries them.** |
| `355015e6` + `19fbff67` + `418f0bb3` + `a8222bd3` | `rad_implicit_ang` (ADI transverse diffusion) was a startup **fatal** on curvilinear meshes | implemented for the cubed sphere: diagonal metric implicit, cross term `g^{23}` explicit (`|cos a| <= 1/2`, unconditionally stable), seam faces as isolated two-cell backward-Euler pairs (`rad_adi_seam_w 0.5`), new `rad_adi_cross_iter`, ADI line partition on an open chain of blocks.  Cartesian **bitwise** (incl. the `box_w8` production config).  cs rate error 0.3 % (l=2) / 1.3 % (l=6) at n=32, conservation 1e-6.  **RKL1 (`sts`) blows up when its substage cap clamps — on the cs use `rad_ang_solver = adi`.**  Test: `tst/test_suite/rad/test_rad_cs_implicit_ang_cpu.py`; design note `docs/dev/cs_implicit_transverse.md`.  Spherical polar keeps its fatal (pole rows). |
| `e4bcdddb`, `0667516d` | three radial-mesh diagnostics in `two_stream`/`red_giant` that were still plane-parallel; three fixes found by the smoke gate | — |

Neither production is affected by any of these: the B and He box runs are plane-parallel Cartesian
and every change was gated bitwise against them.

**Build lesson (hipcc):** dividing by an exact `1.0`, or hoisting an expression into a branch,
changes FMA contraction and breaks bitwise tests.  Append curvilinear forms as *overwrites* of the
verbatim Cartesian statements.

## 4. THE BLOCKER: the spherical two-stream is wrong in the diffusion limit

### 4.1 The term

`1159a8f3` transports the area-weighted intensity against the area-weighted source,
`mu dJ/dr = -kappa rho (J - A B)`.  Substituting `J = A I`:

```
mu dI/dr = -kappa rho (I - B) - 2 mu I / r,
```

i.e. the curvature term `(1-mu^2)/r dI/dmu` has been replaced by a **sink**.  In the transparent
limit that sink gives `I ~ 1/r^2` and conserves `L = 4 pi r^2 F` — exact, and gated by
`tests_r2/rg1d`.  In the **diffusion** limit the same equation gives

```
F_2s = -(mu/(kappa rho)) [ dB/dr + 2 B/r ]   ->   F_2s / F_exact = 1 - H_T/(2r),
                                                  H_T = T/|dT/dr|,  B ~ T^4.
```

The exact first moment has no such term: with Eddington factor `f = 1/3` the sphericity term
`(3f-1)E/r` vanishes identically, and the `r^2` dilution belongs to the **zeroth** moment alone.
So the scheme is O(1) wrong wherever `H_T ~ r`, and returns a flux of the **wrong sign** where
`H_T > 2r`.

### 4.2 The unit test (`tests_r2/thick`, commit `76e21c7d`)

`inputs/tests/two_stream_sph_thick.athinput`: a shell with `rho ~ r^-n`, constant opacity
(`problem/kappa_const`), the **exact radiative-equilibrium** `T(r)` supplied through `ic_profile`,
blend off, one cycle, mode 3 `pcr` — the production operator.  8 cases; `F_2s/F_req` read straight
off the `mlt_dump` face table, with Rosseland `F_raddiff/F_req` as the control.

| set | case | n | tau_tot | r_out/r_in | predicted `1-H_T/2r` | measured | residual |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | n0_t100_r2.0 | 0 | 100 | 2.0 | 0.4613 | 0.4612 | -4e-5 |
| A | n1_t100_r2.0 | 1 | 100 | 2.0 | 0.5308 | 0.5308 | -4e-5 |
| A,B,C | n3_t100_r2.0 | 3 | 100 | 2.0 | 0.6339 | 0.6339 | -2e-5 |
| A | n7_t100_r2.0 | 7 | 100 | 2.0 | 0.7434 | 0.7431 | -3e-4 |
| B | n3_t30_r2.0 | 3 | 30 | 2.0 | 0.6044 | 0.6054 | +1e-3 |
| B | n3_t300_r2.0 | 3 | 300 | 2.0 | 0.6423 | 0.6422 | -1e-4 |
| C | n3_t100_r1.1 | 3 | 100 | 1.1 | 0.9070 | 0.9069 | -7e-5 |
| C | n3_t100_r1.5 | 3 | 100 | 1.5 | 0.7210 | 0.7209 | -3e-5 |

Face by face the ratio tracks `1 - H_T/(2r)` from **0.056** at the bottom of the `n = 0` case (a
flux **20x too small**) to 0.90 at its top.  The error is **independent of optical depth** over a
decade and of shell thickness over `r_out/r_in` = 1.1-2.0 — it is purely geometric.  The Rosseland
control returns `L/(4 pi r^2)` to 0.9997-1.0002, so the deficit is the two-stream's, not the state's.
Residual is 2e-5..7e-4 in the bulk and rises only where `dtau_cell < 0.3` (the diffusion limit itself
failing).  **A thin shell hides this bug; a thick one exposes it.**  Gotcha: `srun` swallows a piped
stdin case list — give `srun < /dev/null`.

On the He star itself the same prediction, formed from the code's own dumped `T`, matches the
measured `F_2s/F_raddiff` to **1-4 % at every face from 0.50 R to 0.98 R**, including the 0.10
minimum at 0.65 R (where `H_T/r = 1.7`) and the recovery to 1 at the surface.  Against tau:
0.68 at tau 300, 0.62 at 100, 0.92 at 50, 0.96 at 20 — **the damage is all at tau 100-300, the FeCZ
interior**.  `ck_nquad 4` and `rt_impl_mixed 0` change nothing; that branch is closed.

### 4.3 The smoke-arm history — five arms, five deaths

Smoke grid, 2 MI300A, `apudev`.  The 1-D gate was abandoned first: a super-Eddington,
Prad-dominated FeCZ has **no hydrostatic radiative-equilibrium state** (15-40 % of L must be
convective), so even with the base flux injected the 1-D column runs away radially at 0.47 turnover.

| arm | configuration | died at | first collapsing cell | `L_out/L` | note |
| --- | --- | --- | --- | --- | --- |
| A | `inner_bc = open` | **0.117** turnover | `i = 5`, **on the inner face** | 0.94 -> 0.69, then 141 | its own boundary runaway: `dmass` **+2.1 %**, ln KE1 e-folds in 27 s (170/turnover), 3.8e4 floors at once |
| B | `wall` + `rt_bottom_flux`, `rad_flux_inner = 1.305278e15` | **0.366** turnover | `r/R = 0.672` | 0.938 -> **0.660**, limps to 0.721 | `efloor` doubling for 0.2 turnover; rho at 0.97 R **-43.6 %**, T -18.2 %; `eos_fail = dfloor = fofc = 0` |
| D | B + MLT sub-grid flux (`mlt_alpha 1.5`) | 0.366 turnover | `r/R = 0.672` | 0.78 | dies within 0.6 s of B |
| BH | 20/300 handover + wall + `rad_flux_inner` | **0.26** turnover | `r/R = 0.633`, the **FeCZ base** | 0.93 -> **0.34**, monotone | `dfloor` 1.15e6 |
| DH | 20/300 + `mlt_alpha 1.5` + vpert 1e-3 | **0.21** turnover | `r/R = 0.662` (the handover radius) | 0.94 -> **0.45** | one-cycle detonation: a cell on both floors given `v1 = -9.4e22` cm/s; `L_rad,cut/L = 8.0e3` |

Convection never starts in any of them: `KE_h/KE_1 = 0.3 %` at death, and `vr_rms` at the kappa peak
oscillates at 0.08-0.39 v_MLT with a ~0.1-turnover period — a coherent radial pulsation, not
overturning.  Nothing is resolution-limited: `eos_fail = 0`, `fofc = 0`, `max c2p it = 0`, the cube
vertices are silent.  Arm A's inner boundary is a mass **source** and an energy **sink** and is
unusable as configured; arm B's inner boundary is innocent.

**The handover experiment and why it was reverted.**  The `t = 0` face budget said `rad_tau_lo/hi =
20/300` closed the 0.50-0.67 R hole (0.10-0.48 -> 0.91-0.99), and on that evidence the values were
made the production default (`6d918992`).  The runs they produced — BH and DH — died **40 % sooner**
than the whole-column arms, so `fe429a58` reverted them.  The input keeps the whole-column blend
`1e5/1e6` because it lives longest, not because it is right.  **Methodological lesson: a cycle-0
budget table is a genuine measurement but not a prediction of survival.**

Two traps in that meter, both documented in `tests_3d/handover/README.md`: (1) `mlt_dump` fires at
the first MLT source call, which under `rt_strang` is **after** the dt/2 two-stream pre-step, so
`F_2s` and `T`/`F_raddiff` in the same row are from different instants; (2) the blend-off
configuration is the **drifted** reference, not the handover runs — proved by the frozen-`cfl`
controls (`H1c/H2c/H3c` agree to four digits at every face although their blends put `w` at 1.0, 1.0,
0.69, 0.0 at the same deep face; `M0b` drifts 8.05e-2 while being bit-identical to `M0`).  On the
frozen state, Rosseland carries L to within 0.7 % over 0.51-0.61 R and 0.3 % over 0.94-0.98 R with a
single **16 % trough at 0.815 R — which is the FeCZ, i.e. the convective flux.**  The earlier "H2
overshoots by 31 % near the photosphere" verdict was entirely the drift and is **retracted**.

**The double count in the MLT closure (suspected, not confirmed).**  The shell-mean deficit closure
forms `D = max(0, F_req - F_cond - F_res - F_2s)`, but `F_cond` already carries the blend weight
(`= w F_raddiff` to four digits) while `F_2s` does **not** carry `(1-w)`.  Under a live blend at
`w ~ 0.94` the face looks nearly fully carried, `D` comes out small, and the un-weighted share is
handed in again at the cut — consistent with DH's cut face passing 8000x its share.  **If the
handover is ever retried, fix this first:** subtract `F_cond + (1-w) F_2s`, not `F_cond + F_2s`.
It is a one-line change and it invalidates arm DH, not the handover as such.

## 5. Recommended path

1. **Rewrite the spherical two-stream as moment equations with a variable Eddington factor.**
   Sphericity in the **zeroth** moment only; the first moment closed with `f(r)` so that `(3f-1)E/r`
   vanishes in the diffusion limit and reproduces the `1/r^2` dilution as `f -> 1`.  This is the only
   option that removes the need for a tau handover at all.
2. **Gate it with both unit tests, before anything else runs**: `tests_r2/rg1d` (transparent limit,
   already the gate for `1159a8f3`) and `tests_r2/thick` (diffusion limit, built today, reports the
   factor directly — one cycle on one GPU).  Both must pass; a fix that is right in one limit and
   wrong in the other is what we already have.
3. **Then re-run the arms** on the unchanged smoke grid, in order B, D, BH, DH, gating on
   `L_rad,out/L` within 5 % of 1 and `efloor = 0` through one turnover.  Only when a survivor clears
   one turnover does arm C (`vpert = 1e-2`, three turnovers, looking for `lnKEh` catching `lnKE1`)
   become meaningful.
4. **Then production**: the grid and cost of §2 (`6 x 320^2 x 192`, ~1 h wall and 8 GPU-h per
   turnover on 8 MI300A), with the block decomposition now free to be 4x4 per panel after
   `94c7165d`.  A 20-turnover relaxation is ~1 day on 8 GPUs; the thermal time (157 turnovers) is
   ~1250 GPU-h.  **No production launch without the user.**

Do **not** attempt the production grid before step 2: it would cost ~100x and reproduce the same
collapse at the same 0.37 turnover.

## 6. Open items

- **`src/bvals/bvals_fc.cpp` carries both seam defects** fixed in `bvals_cc.cpp` by `94c7165d`
  (unconditional monotonicity clamp on an extrapolating stencil; x2x3 edge ghosts unresampled).
  Every cs **MHD** run with more than one MeshBlock per panel is affected, including the dhj cs MHD
  productions.  The corner buffers (slots 48-55) are still a plain copy in both.  Fixing `bvals_fc`
  needs MHD gates.
- **`rad_blend_transverse` under a live blend.**  It is `false` in every current input, which is
  *required* under the whole-column blend (`rad_tau_lo = 1e5` makes `w = 0` on every face, so `true`
  would multiply the x2/x3 conductances by zero and switch horizontal radiative transport off
  entirely).  If a live blend is ever adopted, this setting must be re-derived, not inherited.
- **`rt_profile`'s `T` slot is the pgen's `wtemp` slot, not kelvin** — read only *relative* changes
  from it.  Every table in `tests_3d/` obeys this; anything new must too.
- **`mlt_dump` timing under Strang** (§4.3 trap 1): `F_2s` and `T` in one row are from different
  instants.  Fix the dump point or state the caveat every time.
- **`rt_force_verbose` normalises by `g(r_in)`, not `g(r)`** — cosmetic, was being fixed.
- **Nothing is pushed.**  `he4-presn-global` (36 commits), `cs-implicit-transverse`, `cs-seam-4x4`,
  `he4-rg-plumbing` (merged) and `rt-integration` (`f523f4d4`, `67dd65cd`, `7f9f68a7`) are all local
  only.  `4bcdc855` (the `rt_rad_force` arc-length fix) and `94c7165d` (the seam fix) are the two
  that other branches need.  **Ask the user before pushing.**
- **Disk / inodes.**  `/viper/u2/jinma/ATHENAK` is at **652 GB** of a ~1000 GB quota.  The GPFS
  **inode** quota was hit at ~10:00 today (`git worktree add` failed while a 200 MB `dd` succeeded;
  `mmlsquota` is not installed and `quota` prints nothing).  Diagnose with
  `mkdir .itest; for i in $(seq 1 3000); do : > .itest/f$i || break; done` — it failed at 708.
  Fixed by deleting 76 build dirs in frozen worktrees (~800 inodes each); builds are kept only in
  `wt_merge`, `wt_he4`, `wt_he4_adi`.  A handover memory copy costs ~350 inodes.  Bytes were full
  the same morning: `prod_w7` `rst 00044` was truncated mid-write, and the user approved deleting
  `prod_w7` rst 00000-40 (keep the newest 3), the `prod_w5`/`prod_w6` rst dirs, and `he4_rg_gate`
  (~50 GB).  `prod_w7` is back to 28 `rst` in 56 GB — **prune again.**

## 7. The 36 commits on `he4-presn-global` (`git -C bench/wt_he4 log rt-integration..HEAD`)

| commit | what |
| --- | --- |
| `2e7b30a4` | docs: design note for the implicit transverse radiative operator on the cubed sphere |
| `9e07d03a` | red_giant: mode-3 exact column solve, radiative force and operator splitting |
| `75989503` | red_giant: tau-based top sponge, bottom sponge, and the seed's variable and tau window |
| `a9ee4ec6` | red_giant: `rt_surface` and `rt_profile` dumps, as panel and shell dumps |
| `3a9016f9` | inputs: `he4_presn_cs`, the 4 Msun presupernova He star on the cubed sphere |
| `0667516d` | red_giant/he4_presn_cs: three fixes found by the smoke gate |
| `355015e6` | diffusion: the implicit transverse radiative operator on the cubed sphere |
| `d0ab7bac` | red_giant/he4_presn_cs: free the dumps' host mirrors; cbin `coarsen_factor 2` |
| `1159a8f3` | two_stream_rt: spherical dilution — area-weighted intensities and a volume deposit |
| `2b5a1684` | merge `he4-rg-plumbing` |
| `526791b5` | red_giant: `problem/ic_profile`, an externally supplied 1-D stratification |
| `1512b882` | he4_presn_cs: `ic_profile`, `r_in` at 0.5 R, spherical timescales; 1-D gate artefacts |
| `19fbff67` | diffusion: keep the ADI tridiagonal row's Cartesian statements verbatim, + tests |
| `418f0bb3` | tests_adi: record the `rad_adi_cross_iter` scan |
| `a8222bd3` | diffusion: the ADI line partition on an OPEN CHAIN of MeshBlocks per panel |
| `bbcb06a0` | two_stream mode 3: carry the spherical areas RELATIVE to the column's top face |
| `e4bcdddb` | two_stream/red_giant: three radial-mesh diagnostics that were plane-parallel |
| `244b1e21` | tests_m3: the mode-3 spherical gate artefacts |
| `2c59e35f` | merge `cs-implicit-transverse` |
| `e24f9932` | he4_presn_cs: switch the transverse radiative operator to the implicit cs ADI |
| `4bcdc855` | two_stream_rt: the transverse `Prad grad w` force used the COORDINATE spacing |
| `569c65f8` | tests_3d: the cs-implicit-transverse merge check and the 1-D He gate |
| `e5339b94` | he4_presn_cs: list `problem/rt_bottom_flux` explicitly |
| `94c7165d` | bvals: the cubed-sphere seam resample with more than 2x2 MeshBlocks per panel |
| `192d30f3` | tests_3d/arms: the 3-D he4_presn smoke grid, arms A (open) and B (wall + base flux) |
| `8ac3ee31` | tests_3d/arms: keep the two arms' history files |
| `86aff7f9` | merge `cs-seam-4x4` |
| `f4daa23c` | tests_3d/mlt: the MLT sub-grid flux (arm D) — FAIL, and the term it exposes |
| `e9d3780e` | tests_3d/mlt: keep arm D's history and `rt_profile` |
| `6d918992` | he4_presn: hand the deep interior back to radiative diffusion (`rad_tau_lo/hi = 20/300`) |
| `dd06e962` | tests_m3: the production 1-D He run wedges rather than dies |
| `fe429a58` | he4_presn: keep the whole-column blend as the shipped default; 20/300 is under test |
| `a03d983c` | tests_3d/handover: the tau handover FAILS dynamically — both arms worse than the baseline |
| `872c7208` | tests_3d/handover: the frozen-cfl control REFUTES the drift explanation and clears the IC |
| `76e21c7d` | test: optically thick 1-D spherical unit test of the grey two-stream |
| `75b95308` | tests_3d/handover: control 11760604 closes the meter account (and refutes a second prediction) |

Worktrees: `bench/wt_he4` = `he4-presn-global` (the trunk of today's work),
`bench/wt_he4_adi` = `cs-implicit-transverse` (merged in), `bench/wt_seam4` = `cs-seam-4x4`
(merged in), `bench/wt_he4_rg` = `he4-rg-plumbing` (merged in), `bench/wt_merge` = `rt-integration`
(the production binary).

## 8. What to do first tomorrow

1. **Scope and start the moment-equation rewrite** (§5.1).  Read `tests_r2/thick/README.md` §4 and
   `tests_3d/handover/README.md` §1 first — between them they state the defect, the gate and the
   acceptance criterion.  Run `tests_r2/thick` once on the unmodified code to have the "before"
   table in hand (one 15-min `apudev` job, `sbatch sub_thick.sh`, remember `srun < /dev/null`).
2. **Fix the MLT deficit double count** (§4.3) — one line, and it makes arm DH re-runnable.
3. **Decide with the user what to push**, in particular `4bcdc855` and `94c7165d`, which other
   branches need.
4. **Prune `prod_w7/rst`** (28 files, 0.9 GB each, in a 56 GB directory) and re-run the inode loop
   before creating any new worktree or build.
5. Both productions run unattended to 40 turnovers; the B star reaches it first (~2.7 turnovers/day
   at the current chain rate).  No rolling watches.
