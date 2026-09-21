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
transparent shell (§3c), and the stellar beam was not in flux form at all.

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

**The stellar beam.**  `Qb_g` is added to the gas as a volumetric rate
(`:4288-4296`), so its flux form rests on the kernel's own construction: `transw[cc]` is a
**single running scalar per chain**, updated once per cell (`:3934-3940`), so the value at
face `f` that closes cell `i+1`'s deposit is the identical variable that opens cell `i`'s.
With `ck_spherical` the deposit is `(Phi_beam(f_hi) - Phi_beam(f_lo))/V_i` with
`Phi_beam(f) = A(ie+1) (1-albedo) F_star mu0 T(f)`: the beam is a **parallel pencil
intercepted by the top face**, its cross-section does not shrink with r the way a radial
ray bundle's does, so its flux must **not** be divided by r^2 — but using the *local* face
area with an undiluted flux would make a transparent column absorb `(A_top - A_bot)F` out of
nothing.  Conserving the beam's *power* is the only form that is both.

## 0.2 Per-cell identity, production radial grid (`r_out/r_in` of the ck domain = 1.96, area ratio 3.83)

`max_i |V_i dep_i - (Phi_lo - Phi_hi)| / max|Phi|`, thermal part, from the column dump's
own `A_cell`, `V_cell`, `Src_lw`, `F_lw_net` (`budget.py`):

| | max per cell | integrated over the column |
| --- | --- | --- |
| `ck_spherical = false` | **1.80e-02** | **6.26e-01** |
| `ck_spherical = true` | **1.02e-10** | **1.30e-10** |

1e-10 is round-off *at the solver's own precision* — the chain intensities are `RtF`,
single by default.  The `false` residual is smooth, largest at the top and decaying
downward: it is exactly the missing `-2F/r`.

**Beam, per cell:** exact by construction, as argued in §0.1 — the same `transw` scalar is
the "out" of one cell and the "in" of the next, so there is no second evaluation to
disagree with.  The integrated statement, as a fraction of the power the top face
intercepts, `A(ie+1) F_star mu0`:

| | `sum_i V_i Qb_i / (A_top F_star mu0)` |
| --- | --- |
| `ck_spherical = false` | **0.378** |
| `ck_spherical = true` | **0.991** |

0.991 is `(1-albedo)` times the beam that does not survive to the cut; 0.378 is the old
form depositing as if the column's cross-section were the *local* shell area, i.e.
under-heating the dayside by a factor 2.6 on this grid.

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
| 50 | **-39.1 % / -36.3 %** | **-4.6 % / -2.3 %** |
| 75 | **-40.5 % / -37.8 %** | **-4.2 % / -2.0 %** |

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
| stellar beam | slant path `1/mu0`, deposit `F_star w mu0 (T_hi - T_lo)/dz` | none — and it also did not conserve the intercepted power on a shell |
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
* **The stellar beam is carried as a power, not a flux** — §0.1.

Refused at startup: a Cartesian mesh (the switch is geometry and has no meaning there —
verified by inspection only, since no Cartesian pgen reads it), `rt_ck = false` (verified,
fatals), and `rt_layer_legacy = true`.

Code: `src/utils/two_stream_rt.hpp` — flag + doc note `~:950-1000`; startup guards
`~:1615-1655`; `rt_chain_ck` templated on the flag, the two probe passes `~:3690-3840`, the
face mixing `:3852,:3910,:4046,:4072`, the beam rescaling `:3928-3940`, the `A D` frame
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
