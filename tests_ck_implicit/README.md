# `problem/ck_implicit` — a backward-Euler column solve for the correlated-k two-stream

Code: `src/utils/two_stream_column_ck.hpp` (new; the balance, the Jacobian, the Thomas
solve and every parameter's note), the assembly inside `rt_chain_ck`
(`src/utils/two_stream_rt.hpp`, search `ck_implicit` / `jck`), the driver loop
(`two_stream_rt.hpp`, `picket_fence_two_stream_RT`), the apply hooks
(`:ckimp_`, `skip_de`), and `src/pgen/deep_hot_jupiter_rt.cpp` (plumbing only).
Input: `inputs/tests/dhj_ck_implicit.athinput` (= `dhj_ck_spherical.athinput` plus the
knobs, all at their defaults).  Scripts: `gates.sh`, `relax.sh`.

Default is **off**, and off is the scheme this path is supposed to be running.

---

# 0. THE APPLY KERNEL WAS RUNNING AN INCONSISTENT FORM AT HEAD — READ THIS FIRST

The production binary (27ca5b13) applied the correlated-k source as the **plane-parallel**
divergence, `src = -(Ft - Fb)/dx1`, consistent with its plane-parallel sweep.  Commit
1159a8f3 (he4 work, 09-17) made the *shared* apply kernel area-weighted,
`-(Ft A_t - Fb A_b)/V`, for the **grey** centre-to-centre sweep — which does carry the
`r^-2` dilution in its face conditions — and that form then reached the **correlated-k**
path too, where the sweep does *not* carry it unless `ck_spherical` is set.  Differencing
an undiluted `F` with areas manufactures luminosity at the area ratio: on the production
radial grid that is the +0.63 column-budget violation, the 36-40 % `rt_desum` gap and the
top of the atmosphere collapsing onto `pfloor` that `tests_ck_sph/README.md` reports for
`ck_spherical = false`.  It is a merge regression, not the production scheme.

**Fixed here** (`two_stream_rt.hpp`, the `ckpp_` block at the head of the apply kernel):
when the correlated-k path is active and `ck_spherical` is false, the apply uses areas
`1` and volume `dx1`, i.e. 27ca5b13's arithmetic.  `ck_spherical = true` is untouched, and
so is every grey / picket-fence path (the branch is gated on `ck_on && !cksph_`; with
`rt_ck = false`, `ckpp_` is false and the expression is character-for-character the old
one).

Measured, 20 cycles of the production radial grid, the `rt_desum` sweep-to-gas gap of the
**unchanged semi-implicit apply**, two meshblocks reported separately:

| `ck_spherical` | apply | `rt_desum` gap |
| --- | --- | --- |
| false | HEAD (hybrid, area-weighted) | **-36.3 % / -33.1 %** |
| false | fixed (plane-parallel, = 27ca5b13) | **-8.4 % / -4.4 %** |
| true | area-weighted (correct pairing) | **-8.0 % / -4.2 %** |

The two consistent forms agree with each other to within the spread of the diagnostic, and
the hybrid is four times worse than either.  That is the check asked for: the fixed
off-path reproduces the old scheme's behaviour, the budget closes in the plane-parallel
sense, and nothing collapses onto the floor in these runs.

**Consequence for gate (a):** the HEAD reference binary differs from the fixed tree on a
curvilinear mesh with `ck_spherical = false`, *by design*.  Gate (a) is therefore scoped to
the configurations where nothing was meant to change (below).

All the gates below use the two **consistent** forms — fixed plane-parallel
(`ck_spherical = false`) and `ck_spherical = true`.  The hybrid is not gated; it is only
recorded in the table above.

---

# 1. What was built

`problem/ck_implicit = false | true`.  Per column, the unknown is the cell internal energy
and the balance is

```
R_i(e) = e_i - e*_i - bdt * S_i(e) = 0
```

with `S_i` the **full** correlated-k source the apply block already forms: the
`(1-w)`-weighted two-stream flux divergence, the tau-blend handover, and the direct beam.
`R` is evaluated by **re-running the ordinary ck sweep** at the current iterate — one sweep
per Newton pass — so the fixed point is the exact backward-Euler balance whatever the
Jacobian gets wrong.  At convergence the gas is handed *exactly* `bdt S(e^{n+1})`: no
relaxation factor, no per-cell `(1-e^-x)/x` damping, no `rt_de_max` clip.

**Where the Jacobian came from: rebuilt, not recovered.**  The history does carry a
band-path Jacobian (`0d950ec6`, `515cf64d`, `33fb0604`, `b6413928`: `djd/djup/djuo` and
`jac_g` inside `rt_chain`), but that kernel is the **picket fence** — `gamir1/gamir2`,
`beta`, one grey `B[]` split into two channels — on the **staggered** layers, and the
correlated-k kernel `rt_chain_ck` never had one.  The merge of 09-15/16 moved the band
paths to centre-to-centre layers and the accumulation was dropped; `fd0cf4ff` then made
`rt_implicit_column` fatal under `rt_ck` and states in its message that "the correlated-k /
picket-fence band sweep now uses the centre-to-centre layer branch, which never fills the
mode-3 Jacobian".  What survives today in `rt_chain` is the zeroing and the `djd[NC]`
declarations with nothing writing them.  So the Jacobian here was **rebuilt from the grey
mode-3 template** (`rt_layer_w` in `launch_grey_chain`, `BFaceW`) against the actual
correlated-k centre-to-centre sweep.

**Form.**  At frozen opacity the exchange is linear in the band Planck functions, so
`dS_i/dT_j = sum_{b,g} w_g (dSrc_i/dB_{b,j}) (dB_b/dT)_j`.  Each half layer's source is a
convex combination of the two centre Planck functions (`BFace`, weight `BFaceW`) and of the
two half-layer endpoints (the `dt_l/(dt_l+dt_u)` face interpolation), giving six weights;
`dsrc/dB = e0 dI_in/dB - dem/dB` and `dI_out/dB = (1-e0) dI_in/dB + dem/dB`, with the
carried `dI/dB` truncated to the last cell the ray crossed.  Both sweeps contribute; the
`ck_spherical` face-mixing factor `(1 + beta)` and the up ray's area frame change
`A_below/A_above` are carried exactly (they are multiplicative), the probe passes'
`B`-dependence is not.  `dB_b/dT` is taken by one-sided difference on the same tabulated
Planck-fraction interpolant the sweep reads (`rt_pre_opac`, new `ck_dbdt` array).

Rows: `J_{i,i-1} = -bdt jac0_i/cv_{i-1}`, `J_ii = 1 - bdt jac1_i/cv_i`,
`J_{i,i+1} = -bdt jac2_i/cv_{i+1}`, `cv = e/T` (the same `e ~ T` the semi-implicit apply
assumes; Jacobian-only, so it changes the pass count and not the root).  `jac1 < 0`,
`jac0, jac2 >= 0`: an M-matrix for any `bdt`, Thomas without pivoting.

**The direct beam is in the residual and has no Jacobian term**, because at frozen opacity
it is exactly temperature-independent.  It therefore costs nothing extra and works
identically with `ck_beam_sph` false and true (gated below).  Opacities are frozen over the
step; `problem/ck_impl_refresh_kappa = true` rebuilds them (and hence `tau_ray`) every
pass, at one extra `rt_pre_opac` launch per pass.

**Conduction is NOT folded in.**  Unlike the grey `rt_implicit_column = 3`, which merges
the radial radiative diffusion into `Conduction::ImplicitRadialUpdate`, the radial
conduction below/around `ck_pcut_bar` with the `rad_tau_lo/hi` blend stays the separate
operator it is today.  **The 10 bar hand-over therefore stays split**: the blend's handover
flux enters `S_i` explicitly, lagged by one operator, exactly as with `ck_implicit` off.
Folding it in needs the conduction tridiagonal to accept a correlated-k row and is a
second deliverable.

**Cadence.**  `deep_hot_jupiter_rt` calls the RT from `UserSource`, i.e. once per RK stage:
with `rk2` that is **2 ck sweeps per hydro step today**.  There is no once-per-step switch
on this pgen (`rt_col3_once` is a mode-3 knob and mode 3 is refused here).  `ck_implicit`
does not change the cadence; it makes each stage's source implicit at `ck_impl_maxit`
sweeps instead of one — measured **6.3-7.8 sweeps per stage**, i.e. ~13-16 per hydro step.

Guards (fatal at startup): `ck_implicit` needs `rt_ck`, `rt_use_cons`, `!rt_explicit`,
`!rt_layer_legacy`, `rt_implicit_column == 0`, `rt_outer_iter == 1`.

---

# 2. Gates

Serial CPU, `inputs/tests/dhj_ck_implicit.athinput`: the production radial grid
`x1 = 9.44e9 .. 2.0556e10` (`nx1 = 128`, polynomial stretch, point-mass gravity,
general EOS, radiative conduction + the tau blend), `nx2 = nx3 = 8`, two meshblocks, i.e.
128 columns spanning day, terminator and night.  `build_ckimp_ref` = HEAD unmodified,
`build_ckimp_new` = this tree.

## (a) default off

| configuration | artefacts | result |
| --- | --- | --- |
| `ck_spherical = true`, 20 cycles | `dhj.hydro.hst` + 22 `bin` payloads | **byte-identical** to the reference |
| `ck_spherical = false` (hybrid), 20 cycles | hst + 22 bin + 2 rst payloads | **byte-identical** to the reference *before* the section-0 fix |
| `ck_spherical = false` (fixed) | — | **deliberately different**: see section 0 |
| grey / picket-fence (`rt_ck = false`) | — | by inspection: `ckpp_ = ck_on && !cksph_` is false, the expression is unchanged |

So `ck_implicit = false` is bitwise inert everywhere; the only behaviour change in the tree
is the section-0 correction, which is confined to `rt_ck && !ck_spherical`.
(`rt_ck = false` cannot be run from this input — the pgen fatals, `rad_tau_hi > 0` requires
`rt_ck` — so that row is inspection, not measurement.)

## (b) the sweep-to-gas gap

`rt_desum` (`sum(de dx) / sum(src dt dx) - 1`) for the semi-implicit apply, against the
same ratio formed from the converged implicit state (`ckdesum`, printed by
`problem/ck_impl_verbose`).  20 cycles, cycle 19, two meshblocks separately.

| `ck_spherical` | `ck_beam_sph` | semi-implicit (today) | `ck_implicit` | Newton residual |
| --- | --- | --- | --- | --- |
| false (fixed) | false | -8.38 % / -4.40 % | **+8.6e-09 / -1.8e-06** | 2.9e-08 / 5.4e-09 |
| false (fixed) | true | -8.29 % / -4.35 % | **+1.4e-08 / -2.9e-08** | 7.1e-08 / 5.0e-09 |
| true | false | -8.03 % / -4.20 % | **+1.8e-09 / -1.7e-06** | 2.4e-08 / 5.3e-09 |
| true | true | -8.00 % / -4.18 % | **+2.9e-09 / -1.7e-06** | 2.6e-08 / 4.9e-09 |

All four combinations converge with **0 capped cells and 0 fallbacks**.  The gap is at the
Newton tolerance, which is the target.

**On columns.**  The residual printed is a **max over every column of the rank** — day,
terminator, twilight and night together, including the `ck_beam_sph` twilight columns where
the beam is the only heating.  That is a stronger statement than four named columns would
be, and it is why no per-column table is given: no column is above `7e-08`, so none of them
can be hiding a gap.  The `ckdesum` ratio is the volume-weighted (`dx`) integral over the
same set.

## (c) the column budget in a transient

The section-0.1 result of `tests_ck_sph/README.md` stands unchanged: every face value
`Phi_f = A_f F_f` is written exactly once and read by both neighbours, so
`sum_i V_i S_i = A_bot F_bot - A_top F_top + absorbed beam` *identically*, to the last bit,
for both instantiations.  What was not guaranteed is that the **gas** receives that sum.
The measurement of that is gate (b), and the strongest transient available here is the
first cycles from the analytic initial condition, which is far from radiative equilibrium
(the semi-implicit arm's gap at cycle 0-19 is the -4 to -8 % above; at cycle 0 the implicit
arm needs 20 passes where a settled state needs 4-6).

| | semi-implicit | `ck_implicit` |
| --- | --- | --- |
| `d(sum e dx)/dt` vs `sum S dx`, cycle 19 | off by 4-8 % | off by `< 2e-06` |

**Not done as specified:** the brief asked for a *deliberately perturbed* T profile (+20 %
in a slab).  There is no hook in `deep_hot_jupiter_rt` to perturb the initial temperature,
and adding one is outside the file list for this deliverable, so the cold-start transient
was used instead.  It is a real transient and the numbers are per-cycle, but it is not the
controlled slab perturbation that was asked for.

## (d) the two schemes' relaxed states

`relax.sh`: 600 cycles from the analytic initial condition with `ck_spherical = true`, one
arm semi-implicit and one `ck_implicit`, then three column dumps out of each final restart
(night `mu0 = -0.92`, day `mu0 = +0.92`, day `mu0 = +0.38`), the same pattern
`tests_ck_sph/rerun_dumps.sh` uses.  `max |T_implicit - T_semi|` by pressure range:

| column | `p > 10 bar` | `1-10 bar` | `0.01-1 bar` | `1e-4 .. 1e-2 bar` | `p < 1e-4 bar` |
| --- | --- | --- | --- | --- | --- |
| night, `mu0 = -0.92` | 0.05 K | 0.05 K | 25.3 K | 102.9 K | 251.7 K |
| day, `mu0 = +0.92` | 0.02 K | 0.07 K | 3.3 K | 22.8 K | 111.9 K |
| day, `mu0 = +0.38` | 0.03 K | 0.05 K | 15.8 K | 17.7 K | 138.3 K |

Below 1 bar the two schemes are the same profile to 0.07 K out of 4000-5300 K, i.e. 2e-5
relative: where the cells are thick and the source is not stiff, the semi-implicit apply is
consistent and lands on the same fixed point, as it should.  Above 0.01 bar they separate,
by up to 250 K out of ~1900 K on the night side — which is exactly the thin, stiff region
where the apply's `(1-e^-x)/x` damping and `rt_de_max` are active and where gate (b)'s 4-8 %
of the source is being swallowed.

**Two caveats, and they are not small.**  (i) 600 cycles is **not relaxed** — the two arms
end at `t = 1.4985e4 s` and `1.4923e4 s` (their dt histories differ by 0.4 %), so part of
the upper-atmosphere difference is a time offset and this table cannot separate the two.
(ii) The brief asked which scheme violates its own fixed point if they disagree; that is
*not* answered here.  What can be said is that the implicit arm satisfies its own stated
balance to 5e-9 (gate b) while the semi-implicit arm's `rt_desum` says it misses its own by
4-8 % in the same cells, so the disagreement sits where one arm is known to be
inconsistent.  A proper answer wants a genuinely relaxed pair (several thousand cycles) and
was not run.

## (e) robustness at large dt

Hydro left as it is and `time/cfl_number` raised, so `bdt` scales: **28.98 s** (1x),
289.8 s (10x), 2898 s (100x).  `ck_impl_maxit = 20`, first RT call of the run (the hardest:
the cold start).  `capped` counts cells that hit `ck_impl_dtmax` (per-pass, 0.5) or
`ck_impl_demax` (total excursion from `e^n`, 0.5 — the same bound `rt_de_max` imposes
today).

| `ck_spherical` | dt | passes | residual | capped | fallbacks |
| --- | --- | --- | --- | --- | --- |
| false (fixed) | 1x | 20 (then 4-8 on later calls) | 1.6e-05 | 0 | 0 |
| false (fixed) | 10x | 20 | 8.7e-01 | 8 | 0 |
| false (fixed) | 100x | 20 | 6.5e+01 | 64 | 0 |
| true | 1x | 20 (then 4-8) | 1.3e-05 | 0 | 0 |
| true | 10x | 20 | 1.1e+00 | 56 | 0 |
| true | 100x | 20 | 6.9e+01 | 64 | 0 |

**It does not converge at 10x or 100x, and the caps are not the reason.**  Released
(`ck_impl_demax = 0`, `ck_impl_dtmax = 20`, `maxit = 40`) the iteration *diverges*: at 10x
the residual runs 3.3e+05 -> 2.6e+07, at 100x to 3.7e+12.  The cause is the one the apply
block's own note names ("WHY NOT LINEARIZE ... DIVERGES as E -> 0"): in the optically thin
top cells the emission derivative is ~1e-6 of the source, the linearised row has almost no
diagonal, and Newton asks for `de` many times `e`.  The grey `rt_implicit_column = 3`
handles exactly this with `rt_impl_tau_min` — only cells with `kappa rho dr >~ 1` enter the
tridiagonal, the thin ones keep the bounded nonlinear relaxation.  **That two-level split
is not ported here**, and porting it is the obvious next step.

What the caps do buy: the answer never becomes worse than today's.  `ck_impl_demax` bounds
`|e - e*|` by the same 0.5 `rt_de_max` uses, so a non-converged column is clipped exactly
where the semi-implicit apply would clip it, and the counters say so.  At the production
dt, on both consistent forms, nothing is ever capped.

## (f) cost

100 cycles, same binary for both arms (so this is the cost of the switch alone, not of the
section-0 fix), wall seconds:

| `ck_spherical` | `ck_beam_sph` | off | on | ratio | mean passes/stage |
| --- | --- | --- | --- | --- | --- |
| false (fixed) | false | 23.5 | 119.1 | **5.08x** | 6.35 |
| false (fixed) | true | 25.5 | 137.7 | **5.41x** | 6.63 |
| true | false | 31.3 | 167.5 | **5.35x** | 6.27 |
| true | true | 33.2 | (co-scheduled, not timed) | — | 7.55 |

So a pass costs very nearly one ordinary sweep — the extra `expm1` per half layer and the
two `BFaceW` per face that the assembly adds are in the noise against the `ck_kappa`
look-ups — and the total is the pass count.  **ck sweeps per hydro step: 2 today (one per
rk2 stage), 12-16 with `ck_implicit`.**  Lowering `ck_impl_tol` from 1e-8 to 1e-6 would cut
1-2 passes; that trade was not measured.

---

# 3. Parameters

| knob | default | what it does |
| --- | --- | --- |
| `problem/ck_implicit` | `false` | the switch |
| `ck_impl_tol` | `1e-8` | residual stop, `max_i \|R_i\|/(e_i + eps e_max)` |
| `ck_impl_norm_eps` | `1e-3` | the `eps` above (see `rt_impl_norm` for why) |
| `ck_impl_dtol` | `1e-8` | step stop, `max_i \|de_i\|/e_i`; both must pass |
| `ck_impl_maxit` | `8` | passes = ck sweeps per stage |
| `ck_impl_dtmax` | `0.5` | per-pass cap on `\|de\|/e` |
| `ck_impl_demax` | `0.5` | total cap on `\|e - e*\|/e*`; = today's `rt_de_max` bound |
| `ck_impl_refresh_kappa` | `false` | re-look-up the opacity every pass |
| `ck_impl_verbose` | `false` | the `### ck_implicit` line |
| `ck_impl_debug` | `0` | print this many capped cells with their row |

---

# 4. Files

`gates.sh` (b, e, f), `relax.sh` (d), `a2_*` (gate a, `ck_spherical = true`),
`bpp_*` (gate b, fixed plane-parallel), `b_on_sptrue_*` / `b_off_sptrue_*` (gate b,
spherical), `b_*_spfalse_*` (the HEAD hybrid, section 0), `e2_*` / `e3_*` (gate e),
`f2_*` (gate f), `d_relax_*` / `d_col_*` (gate d).  `.bin` / `.rst` dumps and the build
directories are deleted after measurement (inode quota).
