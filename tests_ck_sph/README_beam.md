# `problem/ck_beam_sph` — the pseudo-spherical direct stellar beam

Default **false** = bitwise the plane-parallel slant path the correlated-k kernel has
always used.  Independent of `problem/ck_spherical` (that one converts the THERMAL
two-stream, this one the BEAM); they compose, and all four combinations are separate
compile-time instantiations of the chain kernel.

## 1. What was wrong

The beam's optical depth to a point was `tau_vertical/mu0` through the column's own
profile, with **`mu0` clamped at 0.1** (`facsw`), and `mu0 <= 0` simply dark.  On a
`r_out/r_in = 2.18` domain that is wrong three ways:

* the secant law diverges at the terminator, and the clamp then silently caps the slant
  depth at 10x the vertical one;
* the true slant depth through curved shells is **finite** at `mu0 = 0` — the Chapman
  function, `Ch -> sqrt(pi X/2)` for `X = r/H` — and is 4.7 (X = 7.6) or 10.3 (X = 61) on
  this atmosphere, not infinity;
* the **twilight** ring `mu0 < 0` is lit on a real planet, down to the radius where the ray
  grazes the opaque body, and received nothing at all.

## 2. The geometry

```
              to star
                 \   theta0                mu0 >= 0 : the ray only climbs.
                  \                        Shells j >= f, weight 1, the target's own
   r_top ----------\------------------     shell contributing from r itself
                    \                      (sqrt(r^2 - b^2) = r|mu0|).
                     \    * target r_f
   ...  ---------------*--\---------------
                       .   \               mu0 < 0 : the ray reaches the target through
                       .    * b  <- tangent  the TANGENT POINT at radius b.  NEAR leg
   r_cut ==============o====\============    r_f -> b, FAR leg b -> r_top, so the shells
   (opaque below; the ray is DARK           between b and r_f are crossed TWICE
    if b <= r_cut = x1f(icut))              (weight 2) and those above r_f once; the
                                            sum starts at the shell containing b.
```

For a target at radius `r` with `cos(theta0) = mu0`, the impact parameter is
`b = r sin(theta0) = r sqrt(1-mu0^2)` and a shell between faces `r_j, r_{j+1}` contributes

```
    ds_j = sqrt(r_{j+1}^2 - b^2) - sqrt(max(r_j^2 - b^2, 0)),
    tau_ray(r) = sum_j w_j ds_j (kappa rho)_j,     w_j = 2 for j < f (mu0 < 0), else 1
```

the `max(.,0)` handling the tangent shell **exactly**.  Above the domain top there is
nothing: `tau_ray(r_top) = 0` (the plane-parallel path instead charges the beam the
unresolved hydrostatic ghost column's depth, which is 2e-4 of the incident flux here).

**Deposit — unchanged in form, local absorption only, no column-power bookkeeping** (see
the retraction in `README.md` §0.1):

```
    Qb_i = (1-albedo) F* w_b (kappa rho)_i [e^-tau(r_i) - e^-tau(r_{i+1})]
                                          / [tau(r_{i+1}) - tau(r_i)]
```

which is **identically** the plane-parallel expression when `tau = tau_vert/mu0` (the
bracket is then `dtau_i/mu0`, and `(kappa rho)_i/(dtau_i/mu0) = mu0/dz_i`), and which never
divides by a chord that disagrees with the optical depths it is paired with.  The
`dtau -> 0` branch uses the same 1e-3 guard the layer coefficients use.

**The flux handed to diagnostics is unchanged**: `Fb_g` is the THERMAL net flux only; the
beam never enters it.  `Q_sw` in the column dump is, as before, the volumetric shortwave
heating rate — but under this switch it is no longer `F* mu0 dT/dz` of a vertical column,
it is the local absorption along the slant ray, and its column integral is **not**
`F* mu0 (1 - e^-tau)` any more (see §3e).

## 3. Implementation

* `mu0` is a pure function of the grid (`rt_pre_geom`: `mu0 = sin(theta) cos(phi)`, the
  substellar point fixed at `phi = 0`); this problem generator has **no non-synchronous
  rotation option**, so `mu0` cannot vary in time and nothing can go stale.  It is
  recomputed every RT call regardless.
* **No chord matrix is stored.**  A per-column triangular `G_ij` would be
  `128*129/2 x 6144 columns x 4 B = 203 MB` on the production cubed-sphere mesh
  (6 x 32 x 32 x 128).  The on-the-fly form costs **one `sqrt` per (face, shell) pair**,
  hoisted out of the g-point loop because `b` does not depend on the chain, and needs no
  memory.  On a GPU that trade is obviously right.
* The one new per-thread array is the `(kappa rho)` column `Krs[NC][NN]`, needed because
  `tau_ray` is not a running sum.  At the production tier (`n1 = 132 -> NN = 136`,
  `NC = RT_NB = 4`, `RtF = float`) that is **2176 B per thread**, on top of the existing
  `I_down` (same size) and, under `ck_spherical`, `Cmx`.  It is sized to **one element**
  when the switch is off — the kernel is templated on the flag, so the default
  instantiation is literally the previous code.
* Extinction shortcut: the face loop breaks once `min_cc tau > 60`.
* GPU-safe: everything the kernel touches is a local or a captured View; `X1F` is the
  existing plane-parallel-aware face-radius helper.

Code: `src/utils/two_stream_rt.hpp` — flag + note `~:1000-1066`, startup guards
`~:1680-1715`, `BSP` template parameter and `Krs` in `launch_ck_chain`, the beam block
after the down-sweep (`// ===== problem/ck_beam_sph: THE PSEUDO-SPHERICAL BEAM =====`),
four-way dispatch at the tier launcher.  `src/pgen/deep_hot_jupiter_rt.cpp` (parameter
only).  `inputs/tests/dhj_ck_spherical.athinput` gains the key.

## 4. Gates

Binaries: `build_cksph_ref` = HEAD (2db7095c) with only the two touched files reverted;
`build_cksph_new` = with them.  Serial CPU, production radial grid (`nx1 = 128`,
`x1 = 9.44e9 .. 2.0556e10`), `nx3 = 16` so more `mu0` values exist, tau blend off.

### (a) default off is BITWISE

20 cycles, hst + bin every cycle + a correlated-k column dump, ref vs new:

| artefact | result |
| --- | --- |
| `dhj.hydro.hst` | **byte-identical** |
| 22 `bin/*.bin` | **data payload byte-identical, 22 of 22** |
| column dump | **bitwise identical** |

### (b) the `mu0 -> 1` limit

`mu0 = 1` is not attainable on this mesh — `mu0` is a grid quantity, and the largest value
the `nx3 = 16` grid offers is **0.9788**.  So the limit is gated two ways:

* **exactly**, in `chapman.py`: at `mu0 = 1`, `b = 0` and the chord formula returns
  `ds_j = r_{j+1} - r_j` with `max|ds_j - dr_j| = 0.000e+00 cm` (`dr = 8.68e7`), so
  `tau_ray` is the vertical depth and the deposit reduces to the old expression term by
  term;
* **in the code**, at `mu0 = 0.9788`: the column-integrated absorbed power is
  `8.5966e+09` (old) against `8.6055e+09` (new), a difference of **+0.10 %** — the pure
  curvature correction at that angle, and nothing else.

### (c) geometry and the Chapman function

`chapman.py`.  **(c1)** The kernel's closed-form chords against an independent reference
that solves `|P + t s| = r_f` for every shell face with a polynomial root finder and
attributes each interval between consecutive crossings to the shell its midpoint falls in
— no `sqrt` identity of the kernel's reused, no marching error.  Over `mu0 = 1, 0.979,
0.435, 0.102, -0.102, -0.195, -0.513` and four target radii:

| | worst relative chord error |
| --- | --- |
| all cases | **1.06e-13** |

Gate was 1e-6.  **(c2)** `tau_ray` for an isothermal exponential atmosphere of constant
opacity on the production shell grid, against the Chapman function computed by its own
quadrature, as the dimensionless enhancement `Ch = tau_ray/tau_vert`:

| `X = r/H` | `theta0` | `Ch` kernel | `Ch` exact | ratio | secant `1/mu0` |
| --- | --- | --- | --- | --- | --- |
| 7.6 (H = 1.24e9, the production H) | 60 | 1.672 | 1.741 | 0.960 | 2.00 |
| | 85 | 3.684 | 3.764 | 0.979 | 11.47 |
| | 90 | 4.676 | 4.749 | 0.985 | **inf** |
| | 95 | 6.139 | 6.221 | 0.987 | **dark** |
| 60.9 (H = 1.55e8) | 60 | 1.966 | 1.923 | 1.022 | 2.00 |
| | 85 | 6.540 | 6.344 | 1.031 | 11.47 |
| | 90 | 10.921 | 10.309 | 1.059 | **inf** |
| | 95 | 19.650 | 20.201 | 0.973 | **dark** |

The residual is the 128-shell discretisation (cell-centre `kappa rho` sampling): 1.3-4 %
at 8 cells per scale height, up to 6 % at 1.8 cells per scale height near the tangent
point, where the path inside one shell is longest.  The point of the table is the last
column: the secant law is **infinite** at 90 deg (the code clamped it at 10) and gives
**nothing** at 95 deg, where the true enhancement is 4.7-10.3 and 6.2-20.2.

### (d) heating profiles, six columns, t = 0

"absorbed per unit top area" is `sum_i Q_i dz_i`, the column integral of the volumetric
rate.  `F* = 8.860e9 erg/s/cm2`, `(1-albedo) = 0.991541`.

| `mu0` | absorbed, old | absorbed, new | `/(F* |mu0|)` old | new | peak-heating `p` old / new [bar] | heating spans `p` [bar] (new) |
| --- | --- | --- | --- | --- | --- | --- |
| **+0.979** | 8.5966e+09 | 8.6055e+09 | 0.9913 | 0.9924 | 1.97e-2 / 1.97e-2 | 4.3e-1 .. 2.1e-7 |
| **+0.435** | 3.8175e+09 | 4.1040e+09 | 0.9911 | 1.0655 | 1.42e-2 / 1.42e-2 | 2.7e-1 .. 2.1e-7 |
| **+0.102** | 8.9502e+08 | 1.6453e+09 | 0.9904 | **1.8206** | 5.57e-3 / 7.60e-3 | 1.5e-1 .. 2.1e-7 |
| **-0.102** | 0 | 5.5673e+08 | 0 | 0.6160 | — / 4.09e-3 | 7.3e-2 .. 2.1e-7 |
| **-0.195** | 0 | 2.5165e+08 | 0 | 0.1459 | — / 2.57e-3 | 3.8e-2 .. 2.1e-7 |
| **-0.513** | 0 | 1.8510e+06 | 0 | 0.0004 | — / 5.69e-5 | 7.3e-4 .. 2.1e-7 |

Read: at `mu0 = 0.979` the two agree to 0.1 %.  By `mu0 = 0.102` the old scheme's clamp has
bitten — it deposits `0.990 F* mu0`, i.e. it still swallows the whole (clamped) beam, but
it puts the peak a factor 1.4 too deep (`5.6e-3` against `7.6e-3 bar`) because the clamped
depth grows too slowly.  **The terminator behaviour is the headline:** the old scheme's
`tau/mu0` diverges as `mu0 -> 0` and is then clamped, and at `mu0 <= 0` it is exactly zero;
the new one lights the twilight ring, with the heating moving steadily **upward** in
pressure as the ray grazes deeper — `4.1e-3 bar` at `mu0 = -0.102`, `2.6e-3` at `-0.195`,
`5.7e-5` at `-0.513` — and dying out smoothly (a factor 300 from `mu0 = -0.10` to `-0.51`)
rather than switching off at a hard edge.

### (e) global absorbed power

Integrating the measured column integrals over `mu0` (`2 pi r^2 d(mu0)`, the six columns
plus `absorbed -> 0` at `mu0 = -1` and the `mu0 = 1` value extrapolated), normalised by
`pi r_abs^2 F*` at `r_abs = 1.2409e10`, the radius where the slant optical depth reaches 1
(measured in `README.md` §0.2):

| | `P / (pi r_abs^2 F*)` |
| --- | --- |
| old (plane-parallel, clamped, dayside only) | **1.0015** |
| new (pseudo-spherical) | **1.1064** |

The old value is 1 to 0.15 %, which is the classical plane-parallel result: a dayside of
opaque vertical columns absorbs exactly `pi r^2 F*`.  The new one is **10.5 % higher**, and
this is worth being explicit about: the pseudo-spherical approximation is **not
photon-conserving**.  Each column is given its own horizontally homogeneous profile, so a
ray that in reality enters through one column's top and is absorbed in a neighbour's
vertical extent is counted in both; near the terminator the per-column integral
`sum Q dz` legitimately exceeds `F* mu0` (1.066 at `mu0 = 0.435`, 1.82 at `mu0 = 0.102`).
That is the known price of the approximation, it is a property of the physics model and not
of this discretisation, and it is why the deposit is defined as a **local absorption rate**
and not by any column-power bookkeeping.  A globally photon-conserving treatment needs
non-radial ray tracing across columns, which is a different (and much more expensive)
scheme.

### (f) cost

100 cycles, production radial grid, `nx2 = nx3 = 8`, serial CPU, whole-run CPU time:

| | `ck_beam_sph = false` | `= true` | ratio |
| --- | --- | --- | --- |
| `ck_spherical = false` | 18.17 s | 20.29 s | **1.117x** |
| `ck_spherical = true` | 25.67 s | 27.75 s | **1.081x** |

Both on against neither: **1.53x**.  The beam's own cost is the `O(N^2/2)` sqrt (one per
face-shell pair, shared by the four chains) plus `O(N^2/2 x NC)` FMA; the extinction break
at `tau > 60` keeps most dayside columns from walking the whole grid.

**Predicted GPU memory: +2176 B of per-thread local storage at the production tier**
(`NC x NN = 4 x 136` floats), no global allocation at all.  For comparison, the rejected
stored-`G_ij` form would have needed 203 MB of global memory and would have had to be
rebuilt whenever `mu0` changed.

## 5. What is not covered

* Exactly `mu0 = 1` is unreachable on this mesh; the limit is gated analytically plus at
  `mu0 = 0.9788` (§4b).
* No constant-opacity hook exists on the correlated-k path (`README.md` §4d), so the
  Chapman comparison (c2) is done on an analytic atmosphere on the same shell grid rather
  than by running the code with a constant `kappa`.
* The slant path assumes the column's own profile on every shell (pseudo-spherical); day/
  night composition and temperature gradients along the ray are not seen.
* Refraction, and the beam's own scattering, are not represented; `albedo` remains a scalar
  attenuation.

## 6. Exact photon budget — gate (e)'s "not photon-conserving" is RETRACTED

`photon_budget.py`.  Gate (e) integrates one 1-D profile used for **every** column, i.e. a
**spherically symmetric** atmosphere.  For such an atmosphere the pseudo-spherical beam is
not an approximation at all: the "column's own profile" *is* the profile, the chords are
the true chords, and the scheme is exact ray tracing.  The double-counting argument in §4e
describes a real limitation of the approximation (day/night **inhomogeneity** is not
seen), but it does not apply to what §4e measures, and it is **not** what gave 1.1064.

### 6.1 What the kernel does with the power that reaches the cut: it is DROPPED

`src/utils/two_stream_rt.hpp:4111` walks faces only down to `icut`
(`for (int i=ie; i>icut-1; --i)`) and the block ends at `:4169` with no term for the
residual transmission `e^-tau(r_cut)`; a `mu0 < 0` ray whose impact parameter falls below
the cut is marked `dark` at `:4121-4122` and `:4150-4152` then deposits exactly zero.
Nothing is deposited at the cut and nothing is handed to the conduction/interior region.
The scheme-consistent reference is therefore

```
    P_exact = F*(1-albedo) int_0^{r_top} 2 pi b a(b) db,
    a(b) = 1 - e^-tau_chord(b)  (b >  r_cut, both legs)
    a(b) = 1 - e^-tau_in(b)     (b <= r_cut, the single incoming leg down to the cut)
```

not `a(b) = 1` below the cut.  On the real profile the difference is 0.003 % (the cut is
opaque), so the choice does not matter here; both are printed.

### 6.2 Profiles

The correlated-k column dumps carry no opacity, so gate (e)'s per-g-point `kappa` is not
recoverable from the artefacts.  Two profiles are used instead.

* **(A) real, effective grey.**  The old plane-parallel deposit is
  `Qb_i = (1-A) F* mu0 (T_{i+1} - T_i)/dz_i`, so the cumulative `sum Q dz` from the top of
  `bd_p979_old/col.txt` — gate (e)'s **own** artefact — gives the band-summed
  vertical transmission `T(r_i)` exactly, hence `tau_vert = -mu0 ln T` and
  `(kappa rho)_j` on every shell.  Real stratification, real total depth (`T` at the
  cut `= 2.04e-4`, `tau_vert(cut) = 8.32`).  It is **grey**: it reproduces the real
  vertical transmission but not the non-grey slowing of decay along long chords
  (quantified in 6.5).
* **(B1/B2) analytic** isothermal exponential on the same 128-shell grid, `X = r/H = 7.6`
  (`H = 1.24e9`, the §4c2 production `H`), `tau_vert(r_abs) = 1`; and a 4-g-point non-grey
  variant (`kappa x 0.03 .. 30`, equal weights) to show the property is **per g-point**.

Both kernels are line-by-line transcriptions: the new beam from `:4102-4169` (chords,
weight 2 below the target, the tangent-shell `max()`, the `dtau > 1e-3` guard, the
`tau > 60` break, the dark rule), the old one from `:4041-4068`.

### 6.3 The budget

`P_scheme = 2 pi int_{-1}^{1} dmu0 sum_i Qb_i(mu0) (r_{i+1}^3 - r_i^3)/3` — the deposit is
per unit **volume**, so the radial weight is `r^2 dr`, not `dz`.  Gauss-Legendre, 400
nodes in `mu0` (converged to 2e-6, table below); `b` quadrature 24-point Gauss-Legendre
per shell.

| profile | `P_new/P_exact` | `P_old/P_exact` |
| --- | --- | --- |
| (A) real, effective grey | **0.999998** | 0.8614 |
| (B1) analytic, `X = 7.6`, grey | **1.000003** | 0.7645 |
| (B2) analytic, non-grey, 4 g-points | **1.000491** | 0.7906 |

The new beam is photon-conserving to **2e-6 (grey) / 5e-4 (non-grey)**.  The residual is
pure discretisation of the deposit formula, not a model property: the pair
`[e^-tau(r_i), e^-tau(r_{i+1})]` belongs to two *different* rays (`b = r sin theta0` is a
function of the face radius), and the flux-difference form divides by the optical-depth
difference of that pair; plus the `dtau <= 1e-3` linearisation.  It is largest for the
non-grey case, where the four chains have very different `dtau` across the same cell.

The **old** beam is not conserving: it loses **14 %** of the power on the real profile
(24 % on the analytic one), the missing photons being the twilight ring it never lights
and the limb chords the `1/mu0` clamp truncates.  It looked like 1.0015 in §4e only
because that normalisation is itself wrong by the same kind of factor (6.5).

Where the photons go, profile (A): **43.5 %** at `b > r_cut` — i.e. on chords that miss
the opaque body entirely, the twilight/limb ring — and 56.5 % at `b <= r_cut`.  Only
0.003 % reaches the cut and is dropped.

### 6.4 The effective absorbing radius — 1.106 (and 1.270) are simply right

```
  r(tau_vertical = 1)        = 1.24003e10 cm        (= README.md's r_abs, 1.2409e10)
  r(tau_chord, limb = 1)     = 1.33762e10 cm        = r(tau_v=1) + 2.65 H
  r_eff = sqrt(P_exact/(pi F*(1-A)))  = 1.39832e10  = r(tau_v=1) + 4.30 H
  (r_eff/r_abs)^2            = 1.2698
```

with `H = 3.68e8 cm` the local pressure scale height at `r_abs` (`X = r/H = 33.7`; the
`H = 1.24e9` of §4c2 is the deep value).  This is the ordinary **transit-radius** effect:
the limb is opaque along grazing chords several scale heights above the vertical `tau = 1`
level, so an extended atmosphere legitimately absorbs more than `pi r_abs^2 F*`.  A number
above 1 is the *expected* answer, and the exact one here is **1.270**, not 1.

### 6.5 So where did 1.1064 come from?

Reproduced exactly by the script (`1.00152` / `1.10639`), and it is a **quadrature and
weighting** artefact, not physics:

* gate (e) puts every absorbed photon at `r_abs`: it integrates `sum Q dz` (per unit top
  area) and multiplies by `2 pi r_abs^2`.  The absorption actually happens over
  `r_abs .. r_top`, so the correct weight `r^2 dr` is larger.  Using the code's own six
  columns with `r^2 dr` instead of `dz r_abs^2` already moves the numbers to
  **old 1.0650, new 1.1825**.
* six `mu0` samples plus a linear extrapolation to `mu0 = 1` under-resolve both the
  dayside peak and the twilight ring.

The remaining gap between 1.1825 and the exact 1.270 is profile (A)'s greyness: a grey gas
with the same *vertical* depth is more opaque along the long twilight chords than the real
correlated-k mixture, where the transparent g-points keep escaping.  The script prints the
comparison column by column — the transcription reproduces the **old** code to
`<= 0.1 %` at all six `mu0` (which validates both the transcription and the reconstructed
`kappa rho`), while for the **new** beam the grey profile absorbs 1.03x (`mu0 = +0.435`),
1.16x (`+0.102`), 1.6x (`-0.102`) and 33x (`-0.513`) the code's non-grey value.  That
spread bounds the real `r_eff` between about 1.18 and 1.27 in units of `pi r_abs^2 F*`;
`P_new/P_exact = 1` is unaffected by it, as (B2) shows directly.

### 6.6 Verdict

`P_new/P_exact = 1.000` on every profile tested.  **The statement in §4e that the
pseudo-spherical beam "is not photon-conserving" is retracted**: for the spherically
symmetric atmosphere that §4e integrates it conserves photons exactly, and 1.1064 was a
low estimate of a number whose true value is ~1.27 for the reason the retraction should
have named — the transit-radius effect.  What remains true, and is worth keeping in §4e in
that form, is the narrower statement: the pseudo-spherical approximation assumes the
column's own profile along the whole ray, so on a *horizontally inhomogeneous* planet the
beam is transported through the wrong medium, and there is then no global photon budget
to appeal to.  That is a statement about day/night contrast, not about the geometry.
