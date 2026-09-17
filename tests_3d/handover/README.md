# The tau handover on the He4 smoke grid: the two-stream's spherical form, and why
# 20/300 is NOT the production default

2026-09-17, viper, branch `he4-presn-global`, binary `../../build_gpu_rg/src/athena`.

`../mlt/README.md` closed on this: at `t = 0` the grey two-stream carries only 0.10-0.74 of
`F_req = L/(4 pi r^2)` over 0.55-0.90 R, while Rosseland diffusion on the **same** state
carries 0.85-1.07, and it recommended putting the optical-depth handover back at a sane
depth. This directory does that, and the answer is: **the diagnosis was right, the fix is
not.**

---

## 1. THE TERM — derived, and confirmed without running anything

Commit `1159a8f3` made the sweep spherical by transporting `J = A I` with source `A B`
(`A` = face area):

    mu dJ/dr = -kappa rho (J - A B).

Substituting `J = A I` gives the intensity equation the scheme *actually* integrates:

    mu dI/dr = -kappa rho (I - B) - 2 mu I / r,

i.e. the radial-ray equation with the curvature term `(1-mu^2)/r dI/dmu` replaced by a sink
`-2 mu I/r`. In the transparent limit that sink reproduces `I ~ 1/r^2`, so `L = 4 pi r^2 F`
is conserved and the scheme is exact -- which is what `../../tests_r2/rg1d` gates.

Two-stream moments (`J+- = A I+-`, `S = (J+ + J-)/2`, `D = (J+ - J-)/2`, `mu = 1/sqrt3`):

    mu dD/dr = -kappa rho (S - A B),      mu dS/dr = -kappa rho D,

so in the diffusion limit `S -> A B` and, since the code divides the reported face flux by
`A` again,

    F_2s = -(mu/(kappa rho)) [ dB/dr + 2 B / r ],

    **F_2s / F_exact = 1 - H_T/(2 r),     H_T = T/|dT/dr|     (B ~ T^4).**

The exact first moment has **no** such term: with the Eddington factor `f = 1/3` the
sphericity term `(3f-1)E/r` vanishes identically, and the `r^2` dilution belongs to the
ZEROTH moment alone. So the scheme is O(1) wrong wherever `H_T ~ r`, and returns a flux of
the **wrong sign** wherever `H_T > 2r`.

**Confirmed on `../mlt/D/mltfaces_D.txt`, no new run needed.** `1 - H_T/(2r)`, formed from
the dumped shell-mean `T(r)` by centred differences, against the measured `F_2s/F_raddiff`:

| r/R | H_T/r | predicted | measured | r/R | H_T/r | predicted | measured |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.506 | 0.555 | 0.722 | 0.686 | 0.861 | 0.857 | 0.572 | 0.548 |
| 0.571 | 1.224 | 0.388 | 0.372 | 0.895 | 0.551 | 0.725 | 0.695 |
| 0.609 | 1.571 | 0.214 | 0.173 | 0.921 | 0.280 | 0.860 | 0.872 |
| 0.648 | 1.729 | 0.135 | 0.100 | 0.942 | 0.121 | 0.940 | 0.961 |
| 0.722 | 1.628 | 0.186 | 0.176 | 0.963 | 0.072 | 0.964 | 0.977 |
| 0.788 | 1.348 | 0.326 | 0.316 | 0.984 | 0.052 | 0.974 | 1.011 |

1-4 % at every face from 0.50 R to 0.98 R, including the 0.10 minimum at 0.65 R and the
recovery to 1 at the surface. **The 2-9x shortfall is this one term** -- not the initial
structure, not the missing convection, and not the column solver.

Against `tau` on this star, `1 - H_T/(2r)` is 0.68 at `tau = 300`, 0.62 at 100, 0.92 at 50,
**0.96 at 20**, 0.97 at 10 and 0.97 at the photosphere. **Above `tau ~ 20` the error is only
3-4 %**; the damage is all at `tau ~ 100-300`, i.e. the FeCZ interior.

## 2. THE t = 0 FACE BUDGET — five one-cycle runs (job 11760033, ~40 s on 2 MI300A)

`budget.sh`. Every run: `time/nlim=1`, all outputs off, `problem/mlt_alpha=1.5` (which is
what makes `mlt_dump` fire), `problem/vpert=0`.
H1/H2/H3 move the handover and let `L` enter at the **conduction wall**
(`inner_bc=wall`, `rt_bottom_flux=false`, `rad_flux_inner=1.305278e15`); M0/Q4 keep the
arm-D whole-column configuration (`rt_bottom_flux=true`, blend at 1e5/1e6) and vary one
two-stream switch each.

**Two bookkeeping facts, both needed to read the dump.** `F_cond` (col 8) already carries
the blend weight -- it is `w*F_raddiff` to four digits at every face. `F_2s` (col 17) does
**not** carry `(1-w)`. So the physical carried fraction is `(F_cond + (1-w) F_2s)/F_req`.
Under a live blend `F_2s = 0` at every face with `w = 1`, because the sweep only runs above
the `icut` face -- correct, not a bug.

### `(F_cond + (1-w) F_2s)/F_req`, every 6th face (R = 2.3717e11)

| r/R | H1 20/300 | H2 5/50 | H3 50/1000 | M0/Q4 blend off | `F_raddiff/F_req` |
| --- | --- | --- | --- | --- | --- |
| 0.506 | 0.993 | 0.993 | 0.926 | 0.735 | 1.071 |
| 0.553 | 0.991 | 0.997 | 0.769 | 0.481 | 1.017 |
| 0.609 | 0.970 | 0.999 | 0.590 | 0.182 | 1.050 |
| 0.667 | 0.910 | 0.963 | 0.476 | 0.107 | 0.987 |
| 0.722 | 0.829 | 0.903 | 0.432 | 0.160 | 0.910 |
| 0.772 | 0.772 | 0.862 | 0.429 | 0.241 | 0.867 |
| 0.815 | 0.742 | 0.844 | 0.457 | 0.336 | 0.851 |
| 0.851 | 0.744 | 0.847 | 0.513 | 0.440 | 0.861 |
| 0.880 | 0.776 | 0.865 | 0.591 | 0.551 | 0.888 |
| 0.902 | 0.838 | 0.906 | 0.698 | 0.681 | 0.929 |
| 0.921 | 0.916 | 0.975 | 0.841 | 0.838 | 0.961 |
| 0.937 | 0.945 | 1.282 | 0.929 | 0.929 | 0.976 |
| 0.952 | 0.956 | 1.089 | 0.957 | 0.957 | 0.986 |
| 0.968 | 0.969 | 0.837 | 0.969 | 0.969 | 0.990 |
| 0.984 | 0.983 | 0.983 | 0.983 | 0.983 | 0.972 |
| 1.001 | 1.185 | 1.185 | 1.185 | 1.185 | 0.797 |

Gate (5 % at all 89 faces with `r/R >= 0.55`): **every run fails.** Worst face: H1 0.739 at
0.828 R; H2 1.309 at 0.942 R (an *excess*); H3 0.425 at 0.757 R; M0/Q4 0.102 at 0.648 R.

* **H1 (20/300) closes the 0.50-0.67 R hole** -- the Schwarzschild-**stable** layer where
  arms A, B and D all put their first collapsing cell -- from 0.10-0.48 to 0.91-0.99.
* The residual over 0.72-0.90 R is **not** a blend artefact: `F_raddiff/F_req` itself dips to
  0.85 at 0.815 R in *every* run, including M0/Q4 where the diffusion operator is inert and
  `F_raddiff` is a pure diagnostic. Those shells are Schwarzschild-**unstable**, so that
  residual is the convective flux.
* **H2 (5/50) overshoots**: it pins `w = 1` out to 0.94 R, into the regime where Rosseland is
  invalid. **H3 (50/1000) is far too deep**: `w = 0.69` already at the bottom face.
* **`rt_impl_mixed = 0` and `ck_nquad = 4` change nothing.** Q4 is bit-identical to the
  default except one round-off digit in one diagnostic column; M0 matches the arm-D baseline.
  The deficit is the spherical form, not the column solver. That branch is closed.

### TWO TRAPS IN THIS METER — read before quoting any number above

1. **The dump is not `t = 0`.** It fires at the first MLT source call
   (`red_giant.cpp:3501-3505`), and under `rt_strang = true` that call happens *after* the
   `dt/2` two-stream pre-step. So `F_2s` (stored by the pre-step) and
   `T_f`/`grad`/`F_raddiff`/`F_cond` (after it) **in the same row are from different
   instants**. Evidence: H1's `F_2s` at `i = 70,76` is bit-identical to M0's while its `T_f`
   there differs by 1.2 %.
2. **Deep `T_f` drifts 0.42 %** between configurations, which -- because `F_raddiff` is a
   gradient -- moves `F_raddiff/F_req` at `i = 4` from M0's 1.0712 to H1's 0.9930, a 7 %
   swing. Leading explanation: with the blend off the radial conduction operator is *inert*,
   with 20/300 it is live, implicit and at `w = 1`, so it relaxes an IC that is over-carrying
   by 7 % toward equilibrium. Either reading brackets 1, so the deep conclusion survives; the
   honest wording is "within 7 % on the IC, within 1 % after one implicit step". Frozen-cfl
   controls (jobs 11760331, 11760604) were queued to settle this.

**Configuration warning.** M0/Q4 took `rad_tau_lo/hi = 1e5/1e6` from the input file. Pass
those explicitly when reproducing the reference; do not rely on the file's defaults.

## 3. THE ARMS — both FAIL, and both worse than the whole-column baseline

`chain.sh`, same smoke grid as `../arms/` and `../mlt/` (nx1 = 96, nx2 = nx3 = 32,
meshblock 96x16x16, 24 blocks, 2 MI300A, `apudev`). Both on the 20/300 handover with
`inner_bc = wall`, `rt_bottom_flux = false`, `rad_flux_inner = 1.305278e15`, `vpert = 1e-3`.
Parameter records were read back from the restart files, not assumed from the command line.

| | arm B (1e5/1e6) | arm D (1e5/1e6 + MLT) | **BH** (20/300) | **DH** (20/300 + MLT) |
| --- | --- | --- | --- | --- |
| `mlt_alpha` | 0 | 1.5 | 0 | 1.5 |
| died at | 1723 s = **0.36** turn | 1724 s = **0.36** turn | 1242 s = **0.26** turn | 998 s = **0.21** turn |
| first collapsing cell | r/R = 0.672 | r/R = 0.672 | r/R = **0.633** (FeCZ base) | r/R = 0.662 |
| `dfloor` | 0 | 0 | **1 149 984** | **627 264** |
| `efloor` | 1230 | 720 | 129 921 | 18 301 |
| `fofc` | 0 | 0 | 0 | **19** |
| `L_out/L` | 0.66, then **recovers**, holds 0.66-0.72 | 0.78 | 0.93 -> **0.34**, monotone | 0.94 -> **0.45**, monotone |
| rho at 0.97 R | -43.6 % | -21.4 % | -25.3 % | **+9.7 %** |
| T at 0.97 R | -18.2 % | -14.7 % | -23.1 % | -15.5 % |

**Ordering, worst to best: DH 0.21 < BH 0.26 < B 0.36 = D 0.36.** Both handover arms are
worse than both whole-column arms, and adding the MLT sub-grid flux to the handover made it
*worse*, which is the opposite of the prediction the handover rested on.

The one thing that behaved as designed: DH is the only arm that **holds the FeCZ-top density
up** (+9.7 % against BH's -25.3 %). The MLT flux does what it is for; it does not save the
run.

DH's death is qualitatively different from every other arm's -- not a dt taper over hundreds
of cycles but a **one-cycle detonation**:

```
### dt COLLAPSE cycle=84 time=997.783 dtold=4.94273 dt=7.22763e-15
    hydro dt is set by cell (m,k,j,i) = (8,18,7,21) gid = 20
    r=1.57015e+11 rho=1e-13 T=1.8571e+12 p=0.1 cs=1.03983e+06 v=(-9.44933e+22,...)
```

`r/R = 0.662`, the handover radius; the cell sits on **both** floors and is then given
`v1 = -9.4e22 cm/s`. The face budget at that moment reads `L_rad,cut/L = 8.0e3`.

### Suspected mechanism — NOT confirmed

The MLT deficit closure and the live blend are filling the same deficit twice. The closure
forms `D = max(0, F_req - F_cond - F_res - F_2s)`, but `F_cond` is the **blended** flux
(`w*F_raddiff`, verified above) while `F_2s` is the **unblended** two-stream flux. Around
`r/R ~ 0.66` in this configuration `w ~ 0.94`, so the face looks nearly fully carried, `D`
comes out small, and the un-weighted share is handed in again at the cut. The cut face
passing 8000x its share is consistent with that. **If the handover is retried, fix this
double count first**: the closure should subtract `F_cond + (1-w) F_2s`, not
`F_cond + F_2s`.

## 4. VERDICT AND WHAT TO DO NEXT

**The defect is real and diagnosed; the handover is not the fix.** `20/300` was briefly made
the production default (`6d918992`) on the strength of the `t = 0` budget alone; `fe429a58`
reverted the values after BH, and this README records DH. The input keeps the whole-column
blend because it is the configuration that lives longest -- not because it is right.

The `t = 0` budget is a genuine and useful measurement and the `mltfaces.txt` files are worth
keeping. What it is not is a prediction of survival: **a cycle-0 table said 20/300 was a
five-fold improvement at 0.55 R, and the run it produced died 40 % sooner.** That is the
methodological lesson of this directory.

In order:

1. **Fix the deficit double count** (section 3) -- it is a one-line change in the closure and
   it invalidates arm DH, not the handover as such.
2. **The moment-equation rewrite**, which is what actually fixes the term in section 1:
   sphericity in the **zeroth** moment only, first moment closed with a **variable Eddington
   factor** `f(r)` so `(3f-1)E/r` vanishes in the diffusion limit and reproduces the `1/r^2`
   dilution as `f -> 1`. Gate it with both unit tests -- `../../tests_r2/rg1d` (transparent,
   exists) and `../../tests_r2/thick` (diffusion). This removes the need for a handover at
   all, and it is the only option that does.
3. Only then re-run arms B, D, BH, DH on this grid.

## Files

* `budget.sh`, `budget_frozen.sh` -- the one-cycle budget runs and their frozen-cfl controls.
* `chain.sh`, `jobs.txt` -- the arm submissions.
* `H1/ H2/ H3/ M0/ Q4/mltfaces.txt` -- the five `t = 0` face budgets of section 2.
* `BH/ DH/` -- `he4.hydro.hst`, `he4.log` (event log), `rt_profile.bin`, the column dump, and
  DH's own `mltfaces.txt`. `.bin`, `.cbin`, `.rst` and `rt_surface.bin` deleted after
  measuring; every number above was read off the full logs first.
