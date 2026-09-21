# `problem/ck_spherical` — the spherical form of the correlated-k two-stream

Commit 8bca3dfa gave the GREY centre-to-centre sweep and the mode-3 column solve the exact
spherical two-stream form and stated, in `src/utils/two_stream_rt.hpp`, that "the
monolithic (rt_split = false) sweep, the correlated-k kernel and the mode-1/2 Jacobian
still carry the PLANE-PARALLEL divergence … r_out/r_in < 1.3 hot-Jupiter domains, per-cent
level".  That does not cover the deep-hot-Jupiter production grid
(`x1 = 9.44e9 .. 2.0556e10 cm`, `r_out/r_in = 2.18`, area ratio 4.75), where the error is
not per-cent but a factor.

`problem/ck_spherical` (default **false** = bitwise today's arithmetic) fixes it.

---

# 0. ENERGY CONSERVATION IN FLUX FORM — the primary criterion

## 0.1 Single-valued face fluxes

**The source the gas actually receives on this path is already a flux divergence, not the
sweep's layer deposit.**  `rt_src_direct` defaults to `false`
(`two_stream_rt.hpp:593`) and `deep_hot_jupiter_rt.cpp` **never reads it**, so on every
correlated-k run `direct_on == false` and `rt_apply` takes

```
two_stream_rt.hpp:4550-4562   Ft += Fb_g(m,b,i+1,k,j);  Fb += Fb_g(m,b,i,k,j);
                              const Real aft_a = AFC(m,k,j,i+1), afb_a = AFC(m,k,j,i);
two_stream_rt.hpp:4585        src = -(Ft*aft_a - Fb*afb_a)/VLA(m,k,j,i);
```

Cell `i` reads `Fb_g(m,b,i+1,k,j)` as its top face and cell `i+1` reads **the same array
element** as its bottom face, multiplied by **the same** `AFC(m,k,j,i+1)`.  So
`Phi_f = A_f F_f` is literally one number for both neighbours, per band block `b` and
summed: what cell `i` loses, cell `i+1` gains, to the last bit.  This is true with the
switch off as well — the apply kernel was always area-correct.  **What was wrong was the
flux itself**: `Fb_g` carried no r^-2 dilution, so `A F` grew like r^2 through a
transparent shell (§3c).  The beam was, and remains, a local volumetric absorption
rate; see the retraction below.

Each face value is written exactly once per band block, on the **lower** side, as
`Phi_f/A_f`:

| face | line | expression |
| --- | --- | --- |
| cut face `icut` | `two_stream_rt.hpp:3975-3988` | `wfc*(u_below - d_below)`, scale `A(icut)/A(icut) = 1` |
| interior face `i+1` | `:4044-4051` | `wfc*(u_below - d_below)*ACC(i)/AFC(i+1)` |
| top face `ie+1` | `:4070-4077` | `wfc*(u_below - d_below)*ACC(ie)/AFC(ie+1)` |

and the up ray's frame change across the face is taken *from* the conservation law,
`I_up = d_above + (u_b - d_b)*A_below/A_above` (`:4049-4051`, `:4075-4077`), rather than
re-derived from a second mixing constant — see §2 and §0.4.

**x1 MeshBlock boundaries: there are none.**  `deep_hot_jupiter_rt.cpp:571` fatals unless
`mesh/nx1 == meshblock/nx1`; the RT sweeps whole radial columns and the mesh is never
decomposed in x1.  So a column, and every face in it, lives inside one MeshBlock and one
rank.

**The `ck_pcut_bar` handover.**  The sweep never touches `i < icut` (`rt_pre_cut`,
`:2216-2232`), and the deep region is carried by radiative conduction.  Both operators
refer to the **same face index and the same area**: the RT writes `Fb_g(...,icut,...)`, a
flux per unit area at face `icut`; `Conduction::AddIsotropicHeatFluxRadiative` writes a
flux per unit area into `flx.x1f(m,IEN,k,j,icut)`
(`src/diffusion/conduction.cpp:1565-1577`, centre-to-centre `dl = x1v(i)-x1v(i-1)` under
its `curv` flag); and the divergence that consumes both is formed by the RK update, which
**does carry the areas on spherical-polar and cubed-sphere meshes**:

```
src/hydro/hydro_update.cpp:61-65   if (use_cubed_sphere || use_spherical_polar)
                                     divf(i) = (flx1(..,i+1)*area1(m,k,j,i+1)
                                              - flx1(..,i  )*area1(m,k,j,i  ))/volume(m,k,j,i);
src/hydro/hydro_update.cpp:91      (Cartesian only) divf(i) = (flx1(..,i+1)-flx1(..,i))/dx1;
src/mhd/mhd_update.cpp:49-52       identical
```

The implicit radial conduction carries the same areas explicitly
(`conduction.cpp:2061` `af = area1_(m,k,j,i)`, `:2095` `vi = vol_(m,k,j,i)`).  With
`rad_blend_use_2s` the RT's own face flux is copied into `pcond->rad_f2s` at the same index
(`two_stream_rt.hpp:4193-4201`) and consumed at that face (`conduction.cpp:1546-1552`,
`:1574`).  **Exception, stated for completeness:** the RKL super-time-stepping path
`rad_sts_all` is Cartesian-only and refuses curvilinear meshes (`conduction.cpp:1795-1798`,
`:2052`).

**The bottom / internal flux: per unit area at WHICH radius?**  On the correlated-k path
`rad_flux_inner` / `rt_bottom_flux` is a *grey*-path knob and is not read by this sweep.
The internal flux enters as an intensity at the **cut face**,
`I_up = B_b(icut) + sigma T_int^4/pi * PlanckFrac_b` when `ck_int_at_cut`
(`two_stream_rt.hpp:3975-3981`) — i.e. a flux per unit area at `r = x1f(icut)`, **not** at
the top and **not** at `x1min`.  Under the spherical form the luminosity it represents is
therefore `A(icut) sigma T_int^4`, and what emerges at the top is diluted by
`(r_cut/r_top)^2` — a factor 3.83 on this grid.  That is the consistent reading of "the
internal flux enters at the bottom of the radiative domain", it is a real change, and it is
**not** silently rescaled.  When the tau blend is active `ck_int_at_cut` is forced false
(`:2145`) because the diffusion operator then carries `F_int` itself.

**The stellar beam — RETRACTION.**  An earlier revision of this switch (and of this file)
carried the beam as a conserved power, `Phi_beam(f) = A(ie+1)(1-albedo) F* mu0 T(f)`, on
the argument that a column must deposit every watt it intercepts at its top face.  **That
is wrong, and it is retracted.**  A parallel pencil is not confined to a spherically
diverging column: the column widens with r, the beam does not, so rays that enter the top
face leave through the column's sides.  The power intercepted at the top face is therefore
*not* what the column absorbs.  That form made the planet absorb `pi r_top^2 F*`
(`r_top = 2.06e10`) instead of `~pi r_abs^2 F*` with `r_abs` the radius where the slant
optical depth reaches 1 — measured below at `1.24e10` — i.e. it over-heated the dayside by
`(r_top/r_abs)^2 ~ 2.7`.  **The claim in an earlier version of this README that the
plane-parallel form "under-heats the dayside by a factor 2.6" was that same error seen from
the wrong side and is withdrawn.**

The beam is therefore **untouched by `ck_spherical`**: under both instantiations the deposit
is the local volumetric absorption rate `kappa rho F* exp(-tau_ray)`, written as a flux
difference across the cell so that the column integral is `F* mu0 (1 - e^-tau_tot)` by
construction, with the pre-existing slant-path approximation `tau_ray = tau/mu0`
(`two_stream_rt.hpp:3925-3952`).  It is bit-for-bit the same arithmetic in both — verified
in §0.2, where every beam number is identical off and on.  `Qb_g` is added to the gas as a
volumetric rate (`:4288-4296`); its single-valuedness rests on `transw[cc]` being one
running scalar per chain, so the face value closing cell `i+1` is the identical variable
that opens cell `i`.

**Scattered / diffuse beam term:** there is none.  The correlated-k shortwave is a pure
direct beam: `Qb_g` goes straight into the gas energy (`rt_apply`, `:4288-4296`) and is
never added as a source to the thermal two-stream.  The only coupling back is through the
gas temperature and hence the Planck function `Bb_g` on the next call.  Scattering (albedo)
enters only as the scalar `(1-albedo)` attenuating the incident flux, Parmentier+2015
(`get_albedo`); no scattered radiation is re-emitted into the diffuse field.

## 0.2 Per-cell identity and the beam, production radial grid (ck domain area ratio 3.83)

`budget2.py`, t = 0, tau blend off so `icut` follows `ck_pcut_bar`.  Two columns,
`mu0 = 0.922` (m0k5) and `mu0 = 0.382` (m1k3).

**THERMAL**, `max_i |V_i dep_i - (Phi_lo - Phi_hi)| / max|Phi|` and the column budget
`[sum V Src + (A_t F_t - A_b F_b)]/max|A F|` — identical for both columns:

| | max per cell | column budget |
| --- | --- | --- |
| `ck_spherical = false` | **1.80e-02** | **+6.26e-01** |
| `ck_spherical = true` | **2.05e-10** | **+2.62e-09** |

1e-10/1e-9 is round-off at the solver's own precision (`RtF`, single by default).  The
`false` residual is smooth, largest at the top and decaying downward: it is the missing
`-2F/r`.

**BEAM** — every number below is **bit-identical off and on**, which is the check that the
correction landed:

| | `mu0 = 0.922` | `mu0 = 0.382` |
| --- | --- | --- |
| absorbed power `sum_i V_i Qb_i` | 1.31810e+29 | 5.66056e+28 |
| `/ (A_top mu0 F*)` | **0.378** | **0.392** |
| `r` at `tau_slant = 1` | 1.2409e+10 (p = 2.0e-2 bar) | 1.2597e+10 (p = 1.2e-2 bar) |
| `/ (A(tau_slant=1) mu0 F*)` | **1.038** | **1.044** |
| beam flux at the top face `/ (F* mu0)` | 0.99133 | 0.99107 |

`(1-albedo) = 0.991541` for this configuration (Parmentier+2015 at `Teff0 = 3082.6 K`,
`g = 942`), so the downward beam flux entering the top face is `(1-albedo) F* mu0` to
**2.1e-4** and **4.7e-4** respectively; the shortfall is exactly the unresolved hydrostatic
ghost column above the domain (`RTTopDtau`), which is intended and is `T_ghost = 0.9998`,
`0.9995`.  The absorbed power is `1.04 x A(tau_slant=1) mu0 F*` and only `0.38 x
A_top mu0 F*` — i.e. the column absorbs what a disc of radius `r_abs` intercepts, not what
crosses the domain top, which is the coordinator's point and the reason the power form was
withdrawn.

## 0.3 Column budget, and what the gas actually receives

Column budget, `[sum_i V_i Src_i + (A_top F_top - A_bot F_bot)] / max|A F|`:

| | `false` | `true` |
| --- | --- | --- |
| `sum_i V_i Src_i` | -1.2163e+29 | -1.8440e+29 |
| `-(A_t F_t - A_b F_b)` | -3.2504e+29 | -1.8440e+29 |
| **BUDGET** | **+6.26e-01** | **+2.6e-09** |

**The sweep-to-gas gap is separately quantified and it is NOT small.**  The existing
`rt_desum` diagnostic (`two_stream_rt.hpp:4429-4434`, printed at `:5290-5300`) compares
`sum(src*dt*dx)`, what the sweep asked for, with `sum(de*dx)`, what the gas received after
the semi-implicit relaxation, the `rt_de_max` limiter and the floors.  It was unreachable
on this path (`rt_cell_report` / `rt_report_every` were never read by
`deep_hot_jupiter_rt.cpp`); two lines of plumbing were added.  100 cycles on the
production radial grid, the two meshblocks reported separately:

| cycle | `ck_spherical = false` | `ck_spherical = true` |
| --- | --- | --- |
| 50 | -39.7 % / -36.8 % | -5.9 % / -3.1 % |
| 75 | **-40.5 % / -37.7 %** | **-5.6 % / -2.9 %** |

(`f2_desum_*`, with the corrected beam and the production conduction + tau blend; the
earlier `f_desum_*` run, without the blend and with the retracted beam, gave -40.5/-37.8 %
and -4.2/-2.0 % — the gap is insensitive to both.)

So the semi-implicit apply currently swallows 2-5 % of the radiative source even in the
corrected run, and 36-40 % in the uncorrected one (whose source is far larger and far
stiffer, because an undiluted thermal flux drives a huge local imbalance).  **This gap is
independent of the geometry and is not fixed by this change.**  It matters as much as the
geometry and should be attacked next: `rt_semi_implicit`, `rt_relax_sub` and `rt_de_max`
are the knobs, and `rt_desum` is now reachable to measure them.

A closed-boundary total-gas-energy versus boundary-luminosity integral was **not** done:
this problem has user radial boundaries, point-mass gravity and an open outer boundary, so
there is no closed budget to form without building one.  `rt_desum` is the honest
substitute and is reported above.

## 0.4 Fluxes evaluated twice — every one, called out

1. **The face mixing constant `c`.**  The down-sweep forms `c = beta(d_above - u_below)`
   from a probe's `u_below`; if the up-sweep re-forms its own `c` from its exact
   `u_below`, the two sides of the face disagree by `beta` times the probe's lag.
   MEASURED with that (natural, and wrong) construction: the per-cell identity degraded to
   3.1e-03 and the column budget to **8.8 %**.  It is therefore **not** re-formed: the
   up-sweep takes `Cmx[cc][f]`, the number the down ray used (`:4046`, `:4072`), and its
   frame change straight from `A_a(u_a - d_a) = A_b(u_b - d_b)`.  The probe's lag then
   survives only as an O(beta) shift in `S`, which nothing conserves.
2. **The centre-to-centre layers.**  Each half layer lies entirely inside one cell and
   each cell owns both of its halves, so `Src_g(i)` is never a two-cell average and the
   face intensity between the halves is computed once and stored (`I_down[cc][f]`).  No
   double evaluation.
3. **The probe passes (P1, P2) do not accumulate anything** — not `Src_g`, not `Fb_g`, not
   `Qb_g`.  They exist only to supply `u_below`.
4. **`Em_g`**, the relaxation rate, is a per-cell quantity (`2 (wfc/mu) kappa rho B`) and
   is not a flux; it is unchanged by this work.
5. **The tau-blend handover term** `d[wF]/dr` in the `direct_on` branch (`:4571-4578`) is
   formed from the raw face fluxes and is therefore also single-valued — but that branch
   is unreachable on the ck path, see §0.1.

---

# 1. What the correlated-k path did before

`problem/rt_split = false` in the production input is a **lie the pgen corrects**:
`src/pgen/deep_hot_jupiter_rt.cpp:489-497` forces `rt_split = true` whenever `rt_ck` is
set ("The correlated-k solver only exists inside the split path"), so the kernel production
runs is `rt_chain_ck` inside the `if (rt_split)` block, not the monolithic `2stream_rt`
picket-fence kernel at the bottom of the file.  The monolithic kernel is grey-only and is
never reached under `rt_ck`.

Inside `rt_chain_ck`, on the centre-to-centre layers (`rt_layer_legacy = false`, default):

| piece | what it did | geometry? |
| --- | --- | --- |
| thermal sweep, per (band, g-point, angle) chain | two independent formal sweeps, down then up, over half layers split at each face | **none** — plain intensities, no area anywhere |
| face flux `Fb_g(f)` | `wfc*(I_up - I_down[f])` | **none** — no dilution, so `L = A F` grew like r^2 |
| deposit `Src_g(i)` | `wfc/dz*(I_in - I_out)` | plane-parallel divergence; `-2F/r` missing |
| stellar beam | slant path `1/mu0`, deposit `F_star w mu0 (T_hi - T_lo)/dz` | none, and **correctly so** — see the retraction in §0.1; `ck_spherical` leaves this untouched |
| internal flux | intensity at the cut face | none |
| top boundary | unresolved hydrostatic column `p/g` above the domain (`RTTopDtau`) | none |
| OLR / photosphere dump | read `Fb_g`, a flux per unit area at its own face | consumer-side only |

# 2. What `ck_spherical = true` does

8bca3dfa's form, unchanged: inside a cell the frame area `A_i = V_i/dx_i` is constant and
the layer solve is the plane-parallel one per unit area against `B`; all the geometry is in
the face conditions, `S` continuous and `A D` continuous.

Specific to this kernel:

* **The probe is two passes.**  The down-sweep's mixing needs `u_below`.  The grey mode-0
  sweep supplies it with a bare plane-parallel up-sweep; that probe **does not dilute at
  all**, so on a 2.18 domain its `u` at the top is wrong by the whole area ratio and
  feeding it into `c` **over**-dilutes the real up ray (MEASURED: `A F` varied by 26 % over
  the top twelve faces).  So P1 is a plane-parallel *down* probe and P2 an *up* probe that
  carries the face mixing against P1's result.  Neither accumulates; both reuse the
  existing `I_down` column plus one new face column `Cmx`, sized to one element when the
  switch is off (the kernel is templated on the flag, so the `SPH = false` instantiation is
  literally the previous code).
* **`A D` continuity is enforced, not inferred** — §0.4 item 1.
* **The stellar beam is untouched** — §0.1, including the retraction of the
  conserved-power form that the first revision of this switch carried.

Refused at startup: a Cartesian mesh (the switch is geometry and has no meaning there —
verified by inspection only, since no Cartesian pgen reads it), `rt_ck = false` (verified,
fatals), and `rt_layer_legacy = true`.

Code: `src/utils/two_stream_rt.hpp` — flag + doc note `~:950-1000`; startup guards
`~:1615-1655`; `rt_chain_ck` templated on the flag, the two probe passes `~:3690-3840`, the
face mixing `:3852,:3910,:4046,:4072`, the `A D` frame
change `:4046-4051,:4070-4077`; three diagnostic columns in `rt_dumpcol` `~:5330-5400`.
`src/pgen/deep_hot_jupiter_rt.cpp:98-100,431-444`.
`inputs/tests/dhj_ck_spherical.athinput` (new).

---

# 3. Accuracy gates

Binaries: `build_cksph_ref` = HEAD with only the touched source files reverted (the tree
was being edited concurrently by other agents, so a plain HEAD build would not have
isolated this change); `build_cksph_new` = the same tree with them.  Serial CPU.

### (a) default off is BITWISE

20 cycles of the existing CPU regression configuration
(`inputs/mhd/deep_hot_jupiter_rt_eos.athinput`, `nx1=64 nx2=8 nx3=8`, meshblock `64/8/4`,
`rt_ck=true`, hst + bin + rst every cycle, plus a column dump):

| artefact | result |
| --- | --- |
| `dhj.mhd.hst` | **byte-identical** |
| 22 `bin/*.bin` + 2 `rst/*.rst` | **data payload byte-identical, 24 of 24** |
| column dump, the 10 columns both binaries write | **bitwise identical** |

Not byte-identical *in full* for one reason: the embedded parameter dump gains
`ck_spherical = 0 # Default value added at run time` (+128 bytes) and the `header offset=`
line that records it.  Everything past `<par_end>` is identical.

### (b) plane-parallel limit

Thin shell, `x1max = 9.44944e9` (`r_out/r_in = 1.001`), stretch off, `ck_pcut_bar = 1e6`
so the whole column is swept, 1 cycle:

| quantity | `max abs(on-off)/max abs(off)` |
| --- | --- |
| `F_lw_net` (face) | **3.9e-06** |
| `Src_lw` (cell) | **3.9e-06** |

Gate was 1e-3.  (`Q_sw` is identically zero there: at p ~ 200 bar the unresolved column
above the domain extinguishes the beam.)

### (c) dilution — the gate today's form FAILS by the area ratio

Over the 40 faces above `p = 1e-4 bar` (`r = 1.514e10 .. 2.056e10`, **area ratio 1.842**):

| | `A F` max/min | `F` max/min |
| --- | --- | --- |
| `ck_spherical = false` | **1.851** | 1.005 |
| `ck_spherical = true` | **1.006** | 1.831 |

The old form holds the flux *per unit area* constant through a transparent shell and
therefore **manufactures luminosity at exactly the area ratio, 1.851 against 1.842**.

### (d) thick limit with a constant opacity — NOT POSSIBLE, and why

There is **no constant-opacity or grey override reachable from the correlated-k path**.
`problem/kappa_const` exists only in `red_giant.cpp:1514` and `cooling_convection.cpp:112`;
`rad_kappa_fac` is applied only on the grey branch (`two_stream_rt.hpp:2263`) and never to
the ck band opacities (`:2352-2356`); `rt_kappa_frozen` is explicitly refused under `rt_ck`
(`:1453-1463`); the per-layer opacity is `ck_kappa(cklk,...) + kc_g` and neither term has a
hook; `ck_table` is mandatory and two self-tests `std::exit` on failure.  So
`inputs/tests/two_stream_sph_thick.athinput` / `tests_r2/thick` cannot be pointed at this
kernel without adding a knob, which was out of scope.

Closest honest equivalents, both done: gate (b) is the diffusion-limit consistency check in
the plane-parallel limit (the layer arithmetic and the deep `dtau >> 1` half layers are
untouched to 4e-6); and the deep half of the gate-(c) column is optically thick (it starts
at 10 bar) and satisfies the per-cell identity to 1e-10 along with the rest.  What is
**not** gated is the absolute accuracy of the diffusion flux against an analytic
`-(4 sigma/3) dT^4/dtau` through the ck path.  Note that the ck kernel never carried the
`J = A I` substitution 8bca3dfa was correcting, so that scheme's `1 - H_T/(2r)` error never
applied to it.

### (3) "agrees with the injected flux from below" — NIGHT column

Production radial grid AND the production physics below the cut: `isotropic_conduction =
radiative`, `rad_implicit_x1`, `rad_cap_ang = 0.5`, `rad_tau_lo/hi = 30/300`,
`rad_kappa_src = table`, `rad_flux_inner = -1` (so `ck_int_at_cut` is forced false and the
diffusion operator carries `F_int`).  **`sigma T_int^4` is imposed per unit area at
`x1min`** — `conduction.cpp:1533-1537`, `if (i == is && fin != 0.0) flx1(m,IEN,k,j,i) +=
fin` — so the injected luminosity is `L_int = A(x1min) sigma T_int^4 = 4.26139e+25`
(`T_int = 537.565 K`).  4000 cycles, `t = 1.114e5 s` (off) / `1.062e5 s` (on); the column
dumped is `mu0 = -0.922`, restarted from the clean-exit restart.

**How steady:** `d(tot-E)/dt` has fallen from `-5.5e32` at `t = 0` to `+2.9e30` (off) /
`-5.5e30` (on), a factor 100-190, and the total energy has moved by -0.27 %.  The
atmosphere is **not** in global radiative equilibrium — it starts from an analytic
hydrostatic state at `Teq = 2500 K` and is still shedding that heat, so every column
radiates ~10^3 `L_int`.  The gate therefore tests the **shape** of `L(r)`, which is what
the geometry controls, not its absolute normalisation.

`L(r)/L_int` over the faces the two-stream owns alone (`w = 0`):

| region | `ck_spherical = false` | `ck_spherical = true` |
| --- | --- | --- |
| `r = 1.307e10 .. 2.056e10` (area ratio **2.47**) | **869.1 -> 2149.2**, a factor **2.473** | **1003.85 -> 1005.04**, flat to **1.2e-3** |
| `r = 1.360e10 .. 2.056e10` (area ratio 2.28) | 940.9 -> 2149.2, factor 2.284 | 1004.93 -> 1005.04, flat to **1.1e-4** |
| `L(top)/L_int` | 2149.2 | 1005.04 |

The `false` column's `L(r)` grows by **exactly the area ratio** (2.473 against 2.470,
2.284 against 2.283) — the flux per unit area is constant instead of the luminosity.  The
`true` column's `L(r)` is flat to `1e-4` once above the last cells that are still
exchanging energy.  ON is flat; it is not at 1, and the reason is the unrelaxed initial
state, not the discretisation: a genuinely star-off, fully relaxed atmosphere would take
far longer than is affordable here and was not run.

Also visible: in the `false` run every cell above `r = 1.52e10` sits on `pfloor`
(`p = 1.000e-09 bar` exactly) with `T` at 557-561 K; in the `true` run the same cells are
at `p = 2e-9 .. 5e-5 bar` and `T ~ 3900-4400 K`.

### (4) "agrees with the beam from above" — DAY columns

Same runs; `A_top F_top,thermal` against `L_int + P_beam`, both normalised by `L_int`:

| column | | `L(top)/L_int` | target `1 + P_beam/L_int` | ratio |
| --- | --- | --- | --- | --- |
| `mu0 = +0.922` | off | 5387.9 | 3122.5 | **1.725** |
| | **on** | 2899.3 | 3024.6 | **0.959** |
| `mu0 = +0.382` | off | 3694.9 | 1317.1 | **2.805** |
| | **on** | 2018.6 | 1254.3 | **1.609** |

`P_beam` is the same physical quantity in both (the beam is untouched); it differs between
off and on only because the relaxed *states* differ (3121.5 vs 3023.6, 1316.1 vs 1253.3).
The `mu0 = 0.922` column closes to **-4.1 %** with the switch on against **+72.5 %** with
it off.  The `mu0 = 0.382` column is still 61 % over: that column is neither relaxed nor
isolated — this is a 3-D run and horizontal advection moves energy between columns, which
a per-column budget cannot see.  Both are large improvements; neither is a converged
equilibrium, and the residual is dominated by the two effects just named plus the 3-6 % of
the source the semi-implicit apply swallows (§0.3).

### (5) T(p) of the relaxed columns, with the corrected beam

From the same restarted dumps, so this is a **real T(p) in kelvin** (the correlated-k
column dump carries the EOS temperature; the earlier version of this gate used `e/rho`
from the binary output because it had no relaxed column dump).  `dT = T_on - T_off`:

| region | night (`mu0 = -0.92`) | day (`mu0 = +0.92`) | day (`mu0 = +0.38`) |
| --- | --- | --- | --- |
| `p > 10 bar` (below the cut) | -0.4 .. +0.3 K | -0.5 .. +0.1 K | -0.5 .. +0.1 K |
| `1 .. 10 bar` | +0.3 .. +5.1 K (mean +1.9) | -1.8 .. +0.1 K | 0.0 .. +1.3 K |
| `1e-3 .. 1 bar` (photosphere) | -440 .. +133 K (mean **-144**) | -197 .. +472 K (mean **+17**) | -259 .. +450 K (mean **+24**) |
| `< 1e-3 bar` (upper atmosphere) | 0 .. +1615 K (mean +879) | +533 .. +3805 K (mean +2862) | +486 .. +3636 K (mean +2758) |

The deep region is untouched (sub-kelvin), as it must be: it is optically thick and below
the ck cut.  The photospheric band moves by a few hundred K either way with a small mean —
the **sign now differs between night and day** (night -144 K, day +17/+24 K), which the
earlier, beam-over-heated version of this gate could not have shown: there the whole
dayside was inflated by the 2.6x beam error.  The upper atmosphere is again a qualitative
difference and not a percentage: with the switch **off** it collapses onto `pfloor`
(`p = 1.000e-09 bar` exactly) and `T` falls to ~560 K; with it **on** the same cells stay
at `p = 1e-9 .. 5e-5 bar` and `T ~ 3900-4400 K`.  That collapse is driven by the *thermal*
error alone now (the beam is identical in both runs), i.e. by an undiluted thermal flux
irradiating the thin upper atmosphere from below.

### (e) SUPERSEDED

The earlier gate (e) (`e_relax_*`, `relaxcmp.py`) was run with the retracted beam and
without the conduction / tau blend.  Its photospheric numbers (+12.9 % mean, +41.9 % max
in `e/rho`) contained the 2.6x beam over-heating and are superseded by (5) above.  The
directories are kept only so the retraction can be checked.

### (e) relaxation on the production radial grid

3000 cycles (`t = 8.5e4 s` off, `8.0e4 s` on — a simulated day, so the photospheric layers
have relaxed radiatively and the deep ones have not), full physics, hydro live.  The binary
output of this problem carries `(dens, velx, vely, velz, eint)` only — under the tabulated
EOS there is **no temperature in it** — so what is compared is the specific internal energy
`e/rho`, a monotone proxy for T at fixed composition, and the density.  Converting to
kelvin would need the EOS inversion (most of `e` here is H2 dissociation and ionisation,
not heat), so **T in K is not reported**; that is the one thing the task asked for that
this session did not produce.

| region | `p(t=0)` | `d(e/rho)` on vs off | `d(rho)` |
| --- | --- | --- | --- |
| deep, below the ck cut | 220 .. 5 bar | **-0.02 %** | +0.02 % |
| just above the cut | 2.6 .. 0.2 bar | -0.1 .. **-2.4 %** | < 1 % |
| photospheric band | 1e-3 .. 1 bar | mean **+12.9 %**, range -3.0 .. **+41.9 %** | up to -11 % |
| upper atmosphere | < 1e-4 bar | **+42x** | see below |

The upper atmosphere is not a percentage difference, it is a **qualitative** one: with
`ck_spherical = false` every cell above `p ~ 1e-4 bar` has collapsed onto the density floor
(`rho = 5.000e-14` exactly, `dfloor`) and `e/rho` has fallen by a factor 40; with the
switch on the same cells sit 2-3 decades above the floor with `e/rho ~ 2e12` and a smooth
profile.  That is the compound effect of the two errors this change fixes — the thin upper
atmosphere was being irradiated from below by an **undiluted** thermal flux while receiving
only 38 % of the stellar power it intercepts — plus the 36-40 % of that source the
semi-implicit apply then swallowed (§0.3).

Emergent-luminosity consistency of the *relaxed* state was not measured: the one-shot
column dump fires at the first RT call, and re-dumping at the end needs a restart, which
was not run.  The `A F` statement at t = 0 is gate (c).

---

# 4. Cost

Gate (e), 3000 cycles, same grid, same machine, serial: `off` 510.2 s CPU, `on` 742.9 s,
**1.46x**.  That is
the two extra probe passes, each paying a `ck_kappa` lookup per cell and chain.  Caching
the opacity in a third `NC x NN` per-thread column would remove most of it at the cost of
the per-thread footprint; not done.

# 5. Files

`budget.py` (per-cell and column identities from a column dump alone), `relaxcmp.py`,
`a_ref/ a_new/` (a), `b_thin_*/` (b), `c_deep_*/` (c), `f_desum_*/` (§0.3),
`e_relax_*/` (e), `athena_ref` + `REF_COMMIT.txt`.  `.bin`/`.rst` dumps are deleted after
measurement (inode quota).
