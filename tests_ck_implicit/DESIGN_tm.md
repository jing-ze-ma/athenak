# DESIGN: `problem/ck_implicit` on the tm sweep (`ck_sweep_form = 1`), production form

Scope: production `ck_spherical = ck_beam_sph = true`, `ck_sweep_form = 1`
(`inputs/production/deep_hot_jupiter_cs_prod4.athinput:286-288`). Design only; line
numbers are HEAD `886bd677`.
Another agent is turning `ckjacp_`/`ckskip_`/`rt_layer_legacy` into compile-time tags,
so the lines near `two_stream_rt.hpp:2155-2210` and `3354-3400` will move.

---

## 0. Bottom line

* **The linearity holds.** At frozen opacity the tm face flux and deposit are exactly
  affine in the band Planck functions, `S = M B + s0`. The tm reflectance `R` does not
  depend on `B` at all. So pass 1's R-recurrence *is* the factorisation (the LU of the
  column), and every later Newton pass only has to re-run the linear forward and backward
  substitution for `Sc` and the rays: no k-table lookup, no `expm1`, no divides, no beam.
* **Estimated cost of a lean pass: about 0.15-0.25 of a full production sweep**, from
  operation counts (section 2). This is not measured. Phase 3 predicted a similar win on a
  fat kernel and measured 0.88 on CPU and a slight loss on GPU. The lean pass must
  therefore be a separate, slim kernel, and its cost must be measured before anything
  else is built on it.
* **Implicit only pays through accuracy, never through dt, in prod4.** dt is on the pure
  hydro CFL: the dt ratio between CFL 0.3 and 0.15 is 2.02, dt is about 5.3 s, and
  `max_eta` does not bind (`tests_tm_prof_growth/README.md` §8). The semi-implicit scheme
  already passed the split-lag gate there (the growth-rate ratio is 1.057, where a lag
  mode would show 2.0). The accuracy defects on record were all measured at a cold start
  with `bdt` = 29 s, 5.5x the production dt, and on the 4-pass form. **Phase T0 below
  re-measures them at the production state before any code is written.**

## 1. Why tm refuses `ck_implicit` today

The refusal is the FATAL at `two_stream_rt.hpp:1729-1735`, reached when `ck_sweep_form`
is not 0. The brief's ~1532 is the general `ck_implicit` guard at `:1530-1541`
(`rt_ck`, `rt_use_cons`, `!rt_layer_legacy`, `rt_implicit_column == 0`), which tm
satisfies. The note at `:1053-1056` gives the reason.

| piece | where | 4-pass-specific? |
| --- | --- | --- |
| **Jacobian assembly** (`jck`) | defined at `:3860`. Used only in the 4-pass body: top half `:4435-4445`, down sweep `:4473-4502`, cut half down `:4559-4565`, `jup` start `:4692`, cut half up `:4707-4716`, up sweep `:4770-4800`, top half up `:4815-4820` | **Yes, the blocker.** The tm block `:3902-4231` `return`s before any of it (`:4231`), so under tm `ckjac_g` would stay as `ck_jac_zero` left it (`:3356-3364`). `J = I`, and the "Newton" step would be plain Picard. The entries differentiate the probe-lagged face mixing: the `(1+beta)` factor at `:4492-4497` is "the whole of its derivative" only because the probe's `u_below` is frozen, and the up ray's `rat_u` at `:4790-4794` is exact only because `I_down`/`Cmx` are frozen (`two_stream_column_ck.hpp:46-50`). tm has no probe, and its face coupling goes through `R`, so the entries must be re-derived. |
| **Frozen operator** (FOP) | tag dispatch `:4849-4884` | **Partly.** The tm code already honours `ckfus`: `Qb_g` is kept at `:3396`, `krof` loads at `:3833`, the beam is skipped at `:4145` and `:4168`, and the top datum is frozen at `:3502-3520`. But `launch_ck_form` passes `std::false_type` for FOP whenever `FRM != 0` (`:4871-4880`, comment `:4849-4852`), so the tm frozen pass is never instantiated. The stored triple is indexed per (cell, chain) half layer, and that indexing fits tm too: both halves share one triple (`:3972-3974`). |
| **colskip** | `:2456`, `:3359`, `:3390`, `:5290`; verdict at `two_stream_column_ck.hpp:569` | No. It is a per-column early return and does not depend on the form. |
| **seed** (`ck_impl_seed`) | `CkImplJacPass` at `two_stream_column_ck.hpp:255`; seed rows at `:640-672` | No. It reads only `cksrc_g`/`ckem_g` from the apply (`:5391-5421`). tm writes `Em_g` with the same weights (`:4141`). |
| **reuse_jac** | `ckjacp_` at `:2167` | No. It gates whatever assembly exists; it just needs one under tm. |
| **thin split** (`ck_impl_arat`) | `:2528-2551` | No, as long as every off-diagonal entry pointing at a cell carries that cell's `dB/dT`. This holds for any form that is linear in `B`, and so it holds for tm. |

**Conclusion: only the assembly and the FOP dispatch are tied to 4-pass. Everything in
the Newton driver and the column solve works for any form.**

## 2. The structure: tm is affine in B at frozen opacity

Read off the code (`:3914-4160`):

* Pass 1: `R <- (R+b)/(1+Rb)`, `R <- T^2 R` (`:3992-4003`). This involves only the
  transmissions `T = 1-e0` and `b = beta`, so it depends on opacity and geometry alone.
  **`dR/dB = 0` exactly.** This corrects `tests_ck_sweep_form/README.md` §4.3, which says
  `dR/dB` needs its own recurrence.
* `Sc <- (1-b) Sc/(1+Rb)` and `Sc <- T(R q + Sc) + p`: linear in the emissions `p, q`,
  with coefficients fixed by `R`, `T` and `b`. `p` and `q` are `cin/cout` combinations of
  `BFace` endpoints. `BFace` is convex in `B` at frozen kappa (`BFaceW`, `:1289-1300`).
  The bottom datum is `B_cut + Iint`: affine, and `Iint` does not depend on T.
* Pass 2: `fsolve` (`:4035-4052`) gives `u_b = (R(1+b) d_a + Sc)/(1+Rb)` and
  `d_b = d_a (1+b)/(1+Rb) - b Sc/(1+Rb)`, which is linear. The top datum
  `(1-e^-dtau) B_ghost` is fixed for the call. `step` is linear in `(I, s_in, s_out)`,
  and the closure `ua - ubf` (`:4134-4140`) is linear too.

So, for each column and chain, `Src = M B + s0` exactly, with `M` a function of
`(kappa rho, geometry)`. Newton therefore needs `B(T)` re-evaluated per cell (a local
operation) and `M` re-applied. Everything that builds `M` can be computed once per call.

**What a lean re-application needs**, per (cell, chain), stored on pass 0:

* `e0`, `cin`, `cout`: 3 values, already cached privately in `Cc0/Cci/Cco` under CCH=2;
* `R` at the face: 1 value;
* the two `BFaceW` weights and `dt_l/dtc`: 3 values.

Per face and column, it also needs `b` and the area ratios, which the existing
`ACC/AFC/BTF` macros already provide. That is **7 Reals per (cell, chain)**: about 2.4 GB
over 88 chains x ~80 ck layers x 6144 columns, or 1.2 GB per rank. Allocate over the ck
layers, not over `n1`. Phase 3's 4-Real store was 1.4 GB (`README_phase3.md:56-60`).
At ~3 TB/s that is about 0.4 ms of streaming per pass per rank, against a sweep of
22.6 ms per rank (13568 ms / 600 calls, `tests_tm_prof_growth/README.md` §7).

**Operation counts per (cell, chain).** The weights are my own assumptions: add/mul/FMA 1,
double divide ~10, `expm1` ~25, `sqrt` ~10, k-table bilinear gather ~20.

| work | full tm pass, production (`CCH=2`, `BSP=1`) | lean pass |
| --- | --- | --- |
| k-table gather (`krof`, pass 1; pass 2 reads `Kpc`) | 1 -> ~20 | 0 |
| `cofs`: `dtau/mu`, `expm1`, `e0/x` | ~45 | 0 (stored) |
| `BFace` x4 (2 divides each) + `sfv` x2 (1 divide) | ~120 | 6 FMA (stored weights) |
| face: tm Moebius + `fsolve` (3 divides) | ~35 | 4 FMA (stored `1/(1+Rb)`, `R(1+b)/(1+Rb)`) |
| layer and `step` x4 arithmetic, `p/q`, `Src/Fb/Em` | ~60 | ~40 |
| **thermal subtotal** | **~280** | **~50 + 7 streamed loads** |
| BSP beam: chord walk O(N^2)/column, `sqrt` shared by 4 chains, early exit at tau > 60 (`:4168-4229`) | ~100-175 (night/twilight walk from `icut`) | **0** (beam frozen) |
| **total** | **~380-450** | **~50-60** |

The lean pass comes to **~0.12-0.15 of a full pass by operation count**. Adding the
per-pass overheads the implicit loop runs anyway (the `B_b(T)` rebuild in `rt_pre_opac`,
the apply at ~3 % of a sweep, which is 433.6/13568 in §7, and `CkImplStep`) gives
**c_pass ~ 0.15-0.25**.

The beam's share of the real sweep is inferred, not measured. tm alone was 1.36x faster
than 4-pass (`tests_ck_sweep_form/README.md` §4.2), yet tm plus beam costs +1.6 % over
4-pass without beam (`tests_tm_prof_growth` §7). That puts the beam at roughly a third of
the production sweep, but it is a cross-binary comparison and **T1 measures it**.

**The warning from phase 3/4.** The chain kernel is clamped at 128 VGPR with ~1.16
waves/SIMD, so it is latency-bound (`bench/bisect_cost/README.md:65-77`). A frozen pass
compiled as an instantiation of the same kernel kept that footprint: 0.88 of a pass on
CPU (`README_phase3.md:124-128`) and no gain on GPU (`README_phase4.md:162-169`). The GPU
arms were also confounded: they ran on `446b7b30`, which contains the runtime-bool FOP
penalty of +27.7 % (bisect_cost r2), while the "plain" arm `a544a761` predates it. **So
the 4.30x of phase 4 is probably overstated by up to ~1.3x.** That is an inference;
re-measure it at HEAD. **Design rule: the lean pass is its own small kernel
(`rt_chain_ck_lin`), not a tag on `rt_chain_ck`.**

**Cost model.** Per RK stage, in units of one production sweep: `X = 1.1` (pass 0,
building and storing M) `+ (n-1) c_pass`. From §7, the sweep is 72 % of the step and the
rest is 28 %, so the step cost relative to the semi-implicit scheme is `0.28 + 0.72 X`.
At n = 5.7 (phase 4 settled state):

| c_pass | 1.0 (T1: tm + assembly only) | 0.6 (T2: frozen beam + kappa) | 0.2 (T3: lean kernel) | 0.15 |
| --- | --- | --- | --- | --- |
| step cost / semi-implicit | ~4.4x | ~3.1x | **~1.75x** | ~1.6x |

At n = 4 with the lean kernel it would be ~1.5x. **Even in the best case, implicit
costs ~1.5-1.8x per step.**

## 3. Per-component requirements (production: spherical thermal + spherical beam + tm)

| component | T-dependence inside a Newton call | Jacobian | linearity (§2) |
| --- | --- | --- | --- |
| **spherical thermal sweep** (`ck_spherical`) | through `B_b(T)` only | Area factors `rat`, `fsc` and the apply's `A/V` (`:5309-5313`) are constant geometric multipliers. They scale entries as in 4-pass (`:4496`, `:4790`). The tm face adds **O(beta) negative entries**: `d_b` depends on `-b Sc/(1+Rb)`, so emission from below lowers the down ray. With `b ~ 0.006` per face (ln 3.446 over ~100 faces) these are small. **Drop negative entries into the residual** to keep the M-matrix, so Thomas needs no pivoting. It becomes quasi-Newton with the same root, the same choice 4-pass makes for the probe (`two_stream_column_ck.hpp:44-54`). | preserved exactly |
| **spherical beam** (`ck_beam_sph`) | **None at frozen opacity.** `Qb` depends on `Krs = kappa rho` (frozen) and on `Fstar`, `ckswf`, the geometry. There is no Planck term. | **No Jacobian term**, as in 4-pass (`two_stream_column_ck.hpp:36-38`). It enters the residual as a constant source. Compute it on pass 0 only: tm already skips it under `ckfus` (`:4145`, `:4168`, and `Qb` is kept at `:3396`). Its only T-dependence is through the opacity lookup, which stays lagged; `ck_impl_refresh_kappa` would re-run the O(N^2) walk every pass. | beam is `s0` |
| **tm recurrence** (`ck_sweep_form = 1`) | through `Sc` and the rays only | A new nearest-neighbour tangent (§4, T1). `dR/dB = 0`, so only `dSc/dB` (pass 1) and `dd/dB` (pass 2) are carried, truncated to the last two cells as `jdn`/`jup` are today. The coefficients come from `Cc0/Cci/Cco`, so the extra `expm1` of `jcof` (`:3870-3880`) is not paid. | exact: `R` is the factor |
| tau-blend handover `w F` (`:5337-5343`) | linear in B (via F) | Dropped, as today (item iii, `two_stream_column_ck.hpp:50`). The conduction hand-over stays split (`README_phase2.md` §3). | affine |

## 4. Phased plan

**T0: decision gate. No code, 1-2 apudev jobs.** Implicit cannot raise dt in prod4
(§0), so it has to fix a measured accuracy defect at the production state.

* **Change:** none. Run semi-implicit tm+sph+beam from the prod4 restart (or
  `cs_mhd_prod3/rst/dhj.00567.rst` with the prod4 `<problem>` block, as
  `tests_tm_prof_growth` did) with `rt_cell_report`/`rt_report_every`, plumbed at
  `deep_hot_jupiter_rt.cpp:628`. Record `rt_desum` per meshblock and the `rt_de_max` clip
  count (`nclip`, `par_reduce_clip4` at `:5287`).
* **Also:** difference `T(p)` above 0.01 bar between the existing `g_cfl03` and
  `g_cfl015` dumps (`tests_tm_prof_growth/g_cfl0*/bin/`) at equal t. If the semi-implicit
  answer moves by only a few K when dt halves, it is already dt-converged.
* **On record (all cold start, `bdt` ~ 29 s, 4-pass):** semi-implicit gap -4.2 to
  -8.4 % (`README.md` §2b) and -2.9 to -5.9 % (`tests_ck_sph/README.md:176-181`), against
  1e-6-level implicit; `|T_impl - T_semi|` 0.13 K below 1 bar, 9-50 K at 1e-4..1e-2 bar
  after 200 cycles (`README_phase2.md` §5), up to 250 K after 600 (`README.md` §2d).
* **Go** only if, at production dt, `|gap| > ~1 %` of the column source **or** the upper
  atmosphere's T moves by more than ~5 K between CFL 0.3 and 0.15. **No-go otherwise.**
  CFL 0.15 costs 2x and implicit ~1.6-1.8x, so implicit wins only where its fix is worth
  more than a ~1.7x smaller dt.

**T1: tm Jacobian (correctness). Compile-time tag `JAC` next to FRM/FOP.**

* **Change:** a `jck` block inside `if constexpr (FRM == 1)`. Pass 1 carries
  `dSc/dB_i`, `dSc/dB_{i+1}` per chain (face: x`(1-b)/(1+Rb)`; layer: `T(R dq/dB + ·) +
  dp/dB`); pass 2 carries `dd/dB` of the last two cells, with `du_b = [R(1+b) dd_a +
  dSc]/(1+Rb)`. Rows go into `ckjac_g` by `atomic_add` as in 4-pass, negative entries
  to the residual. Lift the FATAL at `:1729` for FRM 1 only (sd and `rt_layer_legacy`
  stay refused).
* **Gates:** (a) bitwise off-path (§5); (b) the tridiagonal against a finite-difference
  probe of the full tm sweep on 3 columns (day 0.92, twilight 0.38, night -0.92),
  nearest-neighbour entries within 1e-6 relative; (c) gap `<= 1e-6` max over 600 calls,
  settled restart, production config, `maxit = 8`, `reuse_jac = 1`, `seed = 2`; (d) pass
  count against the 4-pass arm's 5.67 (tm drops no probe, so it could be lower).
* **Expected cost:** ~4.4x semi-implicit, the same as today. This phase is correctness only.

**T2: frozen beam and opacity under tm.**

* **Change:** instantiate FOP=1 for FRM=1 in `launch_ck_form` (`:4871-4880`). The
  `ckfus` guards already exist in the tm body.
* **Gate:** bitwise against T1, whole state, 20 cycles, the way phase 3 gate (b) was
  done (`README_phase3.md:84-102`). Then measure the beam's share of a sweep on apudev
  with rocprof: one pass at `ckfus` against a full pass.
* **Expected:** c_pass ~0.5-0.7, step ~3x.

**T3: the lean linear re-apply kernel. This is the phase that decides the cost.**

* **Change:** a new slim kernel `rt_chain_ck_lin` that reads the §2 store and runs only
  the `Sc` forward and ray backward substitution. Pass 0 is the full T2 kernel and writes
  the store; passes 1+ launch only the lean kernel plus the `B_b` rebuild. **Build it
  first as a stub on frozen inputs and time it before wiring it into the Newton loop.**
* **Gates:** (a) `Src_g`/`Fb_g` bitwise (or <= 1e-14) against a full T2 pass on the same
  `B`; (b) VGPR <= 64 (`tests_ck_kernel_tags/dis.sh`), i.e. >= 2x the occupancy of the
  128-VGPR chain kernel; (c) apudev, 300 cycles, one binary for both arms: **c_pass <=
  0.25 and step <= 1.9x.** If c_pass > 0.4, stop: the operation-count argument has failed
  again, as in phase 3. Memory +1.2 GB per rank.

**T4: large-dt convergence. Only if T0 said go, or `ck_impl_once` is wanted.**

* **Change:** (i) the diagonally scaled norm `|R_i| / (b_ii (e_i + eps e_max))`,
  `b_ii = 1 + 4 bdt E/e` (`README_phase2.md` §2'; replaces
  `two_stream_column_ck.hpp:555`); (ii) a line search with <= 1 backtrack per pass,
  affordable only because with T3 a trial costs one c_pass, not a sweep; (iii) a
  relaxed-start test at pinned `bdt`. Phase 2's restart was re-limited by the dt
  controls, and no fixed-dt knob was found in `driver.cpp`, so a test knob is needed.
* **Gates:** 10x and 100x `bdt` from the settled state converge within `maxit = 20` with
  0 caps and 0 fallbacks. The 1x gap and pass count are unchanged.
* **Cost:** the norm is free; the line search adds <= 1 c_pass per pass.

**T5: production A/B.** Implicit-on-tm (T3 plus `reuse_jac = 1`, `seed = 2`) against
semi-implicit, from the same restart, for 300 cycles plus a >= 2000-cycle accuracy pair:
gap, `T(p)` per pressure band, and floors. Default stays off. The user decides.

**When implicit-on-tm beats semi-implicit.** Both of these must hold:

1. T0 finds a defect at production dt: gap > ~1 % or T not dt-converged.
2. T3 brings the step cost below the cost of the dt reduction that would fix the same
   defect in the semi-implicit scheme, i.e. <~1.7x against CFL 0.15.

A larger dt is not a route in prod4: the RT does not limit dt. It becomes one only if the
hydro CFL stops binding, for example in a future RT-limited configuration.

## 5. Risks and gates

* **Bitwise off-path.** `ck_implicit = false`, `ck_spherical = ck_beam_sph = true`,
  `ck_sweep_form = 1`: 20 cycles, hst plus 22 bin payloads identical to HEAD
  (`paycmp.py`); the same for `ck_sweep_form = 0` with `ck_implicit = true`; FRM=1/JAC=0/
  FOP=0 instructions and scratch unchanged at tier 136 (15372 / 1005,
  `tests_ck_kernel_tags/README.md` §3). JAC must be a compile-time tag: a runtime bool
  in this kernel cost 27.7 %.
* **Sweep-to-gas gap.** `ckdesum` `<= 1e-6` max and median over 600 calls, production
  config, settled restart. The phase-4 settled numbers to match or beat are max 9.69e-7
  and median 3.4e-7 (`README_phase4.md:137`). The tolerance cannot be loosened: the gap
  is ~1e2 times `tol` (`README_phase2.md` §4).
* **GPU cost on apudev.** One binary for all arms, including a semi-implicit arm (this
  cancels the phase-4 lineage confound); `cycles/s = 300 / cpu time`, restart
  `dhj.00567.rst`, prod4 `<problem>` block; rocprof per-kernel ms for the full pass, the
  lean pass and `CkImplStep`.
* **Rank imbalance.** tm+beam already has a 20 % rank spread (§7). colskip and per-column
  pass counts can widen it. Report rank-1 against rank-0 sweep ms.
* **Truncated calls.** Reuse and seed drifted ~1e-2 on cold-start truncated calls
  (`README_phase4.md` §5). Gate on the settled state, and run the missing `maxit = 20`
  cold-start pair (`p4d.sh`).
* **M-matrix.** Log the count of dropped negative entries. If it is ever more than
  `b` x the diagonal, the Thomas stability argument needs a check.
* **Memory and compile time.** +1.2 GB per rank. Up to +3 chain-kernel instantiations
  per tier (FRM1 x JAC x FOP), plus one lean kernel: about +2 min on the HIP TU.
* **Unchanged caveats.** Hand-over still split; opacity lagged over the step; `RT_FP32`
  cannot use the stored triple for J (`:3867`).
