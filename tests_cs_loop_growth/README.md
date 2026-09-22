# The cubed-sphere field loop's L1(B) grows after t ~ 0.2. What is it?

`tests_cs_vertex_cell/README.md` closed its gate (c) with one item open: on the cs MHD
field-loop test (`cs_test` iprob 11, toroidal field seeded from its vector potential) the
per-region `L1(B)` against the exact **static** solution falls until t ~ 0.125 and then
grows in *every* region, vertex-heaviest — at n = 32, VERTEX 2.68e-5 at t = 0.0625,
2.54e-5 at 0.125, 3.31e-5 at 0.25, 6.39e-5 at 0.5, 1.00e-4 at 1.0. The question was
whether that is

* **(i)** the scheme's own O(h^2) error accumulating secularly (a numerical
  dissipation/drift of the loop that any 2nd-order scheme shows),
* **(ii)** a genuine slow instability (error growing *exponentially*, with a rate that does
  not converge under refinement), or
* **(iii)** a spurious drift launched by the cube-vertex hot spot (the 3-5x Maxwell-stress
  residual there driving a flow that then advects the loop).

**The answer is (i).** The exact solution of iprob 11 is STATIC, the error grows because a
**steady** spurious velocity field of amplitude O(h) advects the loop, and the *growth
rate itself converges at second order* (measured 1.88-2.14 over n = 16 -> 32 -> 64 in all
three regions). Nothing grows exponentially: the kinetic energy is constant in time to
1-11 %, `max|v|` to 0.5 %, the magnetic energy to 1e-5, and `div B` stays at 1e-14. At the
FIXED time t = 2.0 the per-region error is a clean second order (1.94/2.02/1.88 and
1.96/1.99/1.97). **(iii) is refuted**: `max|v|` is *smallest* at the cube vertex
(0.93e-3 against 1.17e-3 in the interior at n = 64), and the vertex band is the one region
whose error *saturates* while the interior keeps growing linearly.

**No knob cures it**, and none is needed. `ppmx`/`wenoz` cut the interior growth rate 1.5x;
`cs_wellbalanced_src` 1.1x; `fofc`, `cs_lowbeta_fallback`, `cs_vertex_fill` and
`cs_vertex_fill_cc` change it by nothing at all (they never fire on this problem). The one
switch that matters is `<mhd>/bs_emf`: turning the GS05 upwind corner EMF OFF (the plain
four-face average) turns this into a **real instability** — `L1` 20x larger by t = 2 and
`KE` growing exponentially at dln(KE)/dt = 1.27, x7.8 over the run. The shipped default is
already the cure for the only instability on this test.

Everything below is serial CPU, Release, no MPI, no Slurm.

---

## What the exact solution is

`cs_test` iprob 11 with `bunif = 0`: `B = b0c*(-y, x, 0)` seeded from its vector potential
(`faces_from_potential = true`, so `div B` is at round-off), balanced by
`p = p0 - b0c^2 R^2`. `curl B = 2*b0c*zhat` is uniform, the Lorentz force is
`-2*b0c^2 R Rhat` and the pressure gradient is `+2*b0c^2 R Rhat`, so the state is an exact
**equilibrium at rest**. There is **no background flow and no advection**: the loop is
supposed to be STATIONARY, `v = 0` forever, and the deviation of the face field from its
t = 0 value IS the error at any time. `B.rhat = 0`, so the reflecting radial walls pass no
flux and the problem is closed.

So a growing `L1(B)` is *not* "the exact solution moved and we measured against the initial
one". It is real error. The question is only its character.

## How it was measured

`CSTestResistCheck` prints the per-region face-field `L1` **once**, at the end of a run, so
a time series needed a separate run per time. Instead, a new history function
`CSTestLoopHistory` (`src/pgen/cs_test.cpp`, enrolled for iprob 11 under
`<problem>/user_hist = true`) writes the *same* norms into the `.hst` at the output
cadence — the same region classification by angle (within `<problem>/conv_nband` cells of
`|xi| = pi/4` or `|eta| = pi/4`; VERTEX = both), the same three exact face expressions,
plus `KE` and `max|v|` per region, total `ME`, and `max|div B| dx / |B|`.

It is a diagnostic only; no solver file was touched. **Validated against the finaliser**:
at n = 16, t = 1.0 the `.hst` gives INTERIOR 8.7404e-05, SEAM 1.6763e-04, VERTEX 3.5512e-04
— identical to every digit of `CSTestResistCheck` on the same run and to
`tests_cs_vertex_cell/README.md` section 6c.

```bash
cmake -B build_csloop -D PROBLEM=cs_test -D CMAKE_BUILD_TYPE=Release
make -C build_csloop -j8 && cp build_csloop/src/athena tests_cs_loop_growth/athena_diag
tests_cs_loop_growth/run_res.sh     # n = 16/32/64 to t = 2.0, hst every 0.0625
tests_cs_loop_growth/run_knobs.sh   # the knob sweep at n = 32
tests_cs_loop_growth/run_amp.sh     # the field-amplitude arm at n = 32
python3 tests_cs_loop_growth/growth.py tests_cs_loop_growth/out/*.user.hst
```

`nghost = 3`, `plm`, `hlld`, `rk2`, `cfl = 0.3`, one MeshBlock per panel, `conv_nband`
scaled 2 / 4 / 8 with n so the bands have a fixed physical width.

---

## 1. The resolution scan: the GROWTH RATE is second order

`L1(B)` per region, normalised by `max|B_exact|`, against time (`plm`, baseline):

| n | region | t=0 | 0.125 | 0.25 | 0.5 | 1.0 | 1.5 | 2.0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 16 | INTERIOR | 4.911e-05 | 6.235e-05 | 5.594e-05 | 5.861e-05 | 8.760e-05 | 1.399e-04 | 1.849e-04 |
| 16 | SEAM | 7.022e-05 | 7.979e-05 | 7.653e-05 | 9.734e-05 | 1.680e-04 | 2.495e-04 | 3.184e-04 |
| 16 | **VERTEX** | 9.894e-05 | 8.739e-05 | 1.008e-04 | 1.920e-04 | 3.557e-04 | 4.268e-04 | 4.612e-04 |
| 32 | INTERIOR | 1.220e-05 | 1.463e-05 | 1.372e-05 | 1.344e-05 | 2.108e-05 | 3.582e-05 | 4.814e-05 |
| 32 | SEAM | 1.758e-05 | 1.929e-05 | 1.979e-05 | 2.561e-05 | 4.314e-05 | 6.204e-05 | 7.869e-05 |
| 32 | **VERTEX** | 2.544e-05 | 2.546e-05 | 3.353e-05 | 6.403e-05 | 1.003e-04 | 1.151e-04 | 1.254e-04 |
| 64 | INTERIOR | 3.038e-06 | 3.613e-06 | 3.440e-06 | 3.346e-06 | 5.317e-06 | 9.005e-06 | 1.234e-05 |
| 64 | SEAM | 4.398e-06 | 4.931e-06 | 5.248e-06 | 6.823e-06 | 1.125e-05 | 1.565e-05 | 1.975e-05 |
| 64 | **VERTEX** | 6.436e-06 | 7.533e-06 | 1.056e-05 | 1.803e-05 | 2.584e-05 | 2.944e-05 | 3.195e-05 |

**Order at the FIXED time t = 2.0** (which is what a convergence statement means):

| region | n=16 | n=32 | order | n=64 | order |
| --- | --- | --- | --- | --- | --- |
| INTERIOR | 1.849e-04 | 4.814e-05 | 1.94 | 1.234e-05 | 1.96 |
| SEAM | 3.184e-04 | 7.869e-05 | 2.02 | 1.975e-05 | 1.99 |
| **VERTEX** | 4.612e-04 | 1.254e-04 | **1.88** | 3.195e-05 | **1.97** |

Clean second order in every region, VERTEX/INTERIOR 2.49 -> 2.60 -> 2.59 (saturating, as in
`tests_cs_vertex_cell` section 6c at t = 1.0).

**The growth rate itself, `dL1/dt` fitted linearly over t = 0.25..2.0:**

| region | n=16 | n=32 | order | n=64 | order |
| --- | --- | --- | --- | --- | --- |
| INTERIOR | 8.152e-05 | 2.215e-05 | **1.88** | 5.745e-06 | **1.95** |
| SEAM | 1.483e-04 | 3.548e-05 | **2.06** | 8.620e-06 | **2.04** |
| **VERTEX** | 2.027e-04 | 4.652e-05 | **2.12** | 1.058e-05 | **2.14** |

This is the decisive table. **A dissipation / truncation drift has a rate that scales as
h^2 for a 2nd-order scheme; an instability does not.** The rate is h^2 to within 7 % in
every region, at both refinement steps.

**It is not exponential.** The `dln L1/dt` fits over the same window (0.78 / 0.85 / 0.73 at
n = 16, 0.86 / 0.80 / 0.58 at n = 32, 0.89 / 0.76 / 0.48 at n = 64) are not a growth rate
at all: the curves are concave, the interior rising close to LINEARLY in t and the vertex
SATURATING. At n = 32 the vertex increments per 0.25 of time run 3.0e-5 (t 0.25->0.5),
0.9e-5 (1.0->1.25), 0.54e-5 (1.75->2.0) — decelerating by 6x. An exponential instability
would give a *constant* `dln L1/dt` that does not fall with n; this one falls with n in the
vertex band (0.73 -> 0.58 -> 0.48) because the phase the run samples moves with n.

## 2. Nothing is growing: KE, ME, max|v|, div B

| n | KE(0.125) | KE(2.0) | dln KE/dt | max\|v\| INT (0.25 -> 2.0) | max\|v\| VTX (0.25 -> 2.0) | dME/ME | max\|divB\|dx/\|B\| |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 16 | 1.971e-05 | 2.190e-05 | 0.051 | 4.689e-03 -> 4.688e-03 | 3.564e-03 -> 3.497e-03 | +8.1e-05 | 5.9e-15 |
| 32 | 2.435e-06 | 2.540e-06 | 0.027 | 2.337e-03 -> 2.337e-03 | 1.836e-03 -> 1.820e-03 | +3.0e-05 | 7.8e-15 |
| 64 | 2.998e-07 | 3.079e-07 | 0.014 | 1.172e-03 -> 1.173e-03 | 9.318e-04 -> 9.283e-04 | +8.6e-06 | 1.3e-14 |

* A spurious velocity field is established in the first few steps and then **stays put**:
  `max|v|` is constant to 0.5 % over t = 0.25..2.0 at every resolution, and the total `KE`
  rises 11 % / 4 % / 3 % over the whole run — and that residual drift *shrinks* with n.
* Its amplitude is **first order in Linf and second order in the L1 sense**: `max|v|`
  4.688e-3 / 2.337e-3 / 1.173e-3 is exactly a factor 2 per refinement (v ~ h), while
  `KE ~ Int rho v^2 dV` falls 8.09x / 8.12x (~ h^3).
* The **magnetic energy is conserved to 1e-5** and *rises* slightly, so the loop is not
  being resistively dissipated: what grows is a *redistribution* of `B`, i.e. the steady
  flow shearing the loop, exactly the `v.grad B` a truncation velocity produces.
* `div B` never leaves round-off (CT is doing its job).

Put together: `L1(B)` grows because `B` is being advected at ~`v_trunc` by a steady
O(h) flow that never changes, so the displacement, and with it the error, accumulates
linearly in t while `|grad B|` lasts. That is definition (i).

## 3. Is the cube vertex driving it? No

| n | region | max\|v\| at t=2 | KE at t=2 | dL1/dt |
| --- | --- | --- | --- | --- |
| 32 | INTERIOR | 2.337e-03 | 1.604e-06 | 2.215e-05 |
| 32 | SEAM | 2.344e-03 | 8.093e-07 | 3.548e-05 |
| 32 | **VERTEX** | **1.820e-03** | 1.266e-07 | 4.652e-05 |
| 64 | INTERIOR | 1.173e-03 | — | 5.745e-06 |
| 64 | **VERTEX** | **9.283e-04** | — | 1.058e-05 |

The spurious flow is **weakest** at the cube vertex — 22 % below the interior at n = 32 and
21 % at n = 64 — and it does not grow there any faster than anywhere else. There is no
vertex-launched flow. The vertex's larger `dL1/dt` (2.1x the interior) is the same 2-5x
Maxwell-stress *constant* `tests_cs_vertex_cell` measured on the static gates, now acting on
the same steady drift; and it is the region whose error saturates first, because the loop
there is displaced only until the local error field and the displacement balance.

The `BY DISTANCE from the panel edge` profile at n = 32, t = 2.0 (d0 1.033e-04, d1
8.836e-05, d2 7.660e-05, d3 6.638e-05, d4 5.811e-05, d5+ 4.624e-05) decays smoothly by only
2.2x over five cells — a broad, bulk error, not a boundary layer at the seam or the vertex.
By radial index it rises monotonically outward (2.87e-05 at i0 to ~5e-05 at the top), which
is where `|B| = b0c r` and `|grad B|` are largest. Both say the same thing: the growth lives
in the BULK, and the vertex only scales it.

## 4. The knob sweep, n = 32, t = 2.0

`dL1/dt` over t = 0.25..2.0, and the ratio to the baseline (lower is better):

| arm | INT | ratio | SEAM | ratio | VERTEX | ratio | dln KE/dt |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **baseline** (plm, all defaults) | 2.215e-05 | 1.00 | 3.548e-05 | 1.00 | 4.652e-05 | 1.00 | 0.027 |
| `mhd/fofc = true` | 2.215e-05 | 1.000 | 3.548e-05 | 1.000 | 4.652e-05 | 1.000 | 0.027 |
| `mhd/cs_lowbeta_fallback = 0` | 2.215e-05 | 1.000 | 3.548e-05 | 1.000 | 4.652e-05 | 1.000 | 0.027 |
| `mesh/cs_vertex_fill = false` | 2.216e-05 | 1.000 | 3.550e-05 | 1.001 | 4.696e-05 | 1.009 | 0.027 |
| `mesh/cs_vertex_fill_cc = true` | 2.217e-05 | 1.001 | 3.550e-05 | 1.001 | 4.682e-05 | 1.006 | 0.027 |
| `mhd/cs_wellbalanced_src = true` | 2.000e-05 | 0.903 | 2.987e-05 | 0.842 | 4.084e-05 | 0.878 | 0.019 |
| `mhd/reconstruct = ppmx` | **1.436e-05** | **0.648** | 3.435e-05 | 0.968 | 3.642e-05 | 0.783 | 0.007 |
| `mhd/reconstruct = wenoz` | **1.504e-05** | **0.679** | 3.444e-05 | 0.971 | 3.732e-05 | 0.802 | 0.008 |
| **`mhd/bs_emf = true`** | **4.457e-04** | **20.1** | **4.346e-04** | **12.2** | **4.704e-04** | **10.1** | **1.269** |

* **`fofc` and `cs_lowbeta_fallback` never fire on this problem** (beta ~ 20 everywhere, no
  shocks): identical to the baseline to every printed digit. They are not the cause and
  cannot be the cure.
* **The cube-vertex halo fills are irrelevant here** — turning the face-centred fill OFF, or
  the cell-centred one ON, moves the vertex rate by 0.6-0.9 %. One more confirmation that
  this is a cell-update effect, not a halo effect (same verdict as the earlier "oracle").
* `cs_wellbalanced_src` buys a uniform ~12 % and, unlike on `strat`, does not hurt the
  vertex. Small, real, cheap.
* **`ppmx`/`wenoz` are again the only real lever**, and only in the panel INTERIOR (1.5x)
  and at the vertex (1.25x); they do essentially nothing at the seam. They also cut the
  residual KE drift 4x (dln KE/dt 0.027 -> 0.007), i.e. they make the spurious flow more
  nearly steady. Consistent with `tests_cs_regions` section 3: on `loop` no reconstruction
  improved the vertex L1 at a *fixed short* time — the benefit here is on the *rate*.
* **`bs_emf = true` is the one setting that turns this into a genuine instability**, and it
  is the diagnostic switch, not a default. Replacing the GS05 upwind corner EMF by the
  plain four-face average makes `KE` grow *exponentially* (dln KE/dt = 1.27, KE x7.8 from
  t = 0.125 to 2.0, `max|v|` rising 2.3e-3 -> 5.8e-3 and now LARGEST at the vertex,
  5.80e-3) and `L1(B)` 10-20x worse. This is what (ii) actually looks like, and it is
  absent from the shipped configuration. Keep `bs_emf` off.

## 5. The field-amplitude arm: the drift is driven by the Maxwell-stress flow

n = 32, `b0c` varied (`L1` is normalised by `max|B_exact| ~ b0c`, so these are *relative*
errors):

| b0c | max\|v\| INT | ratio to b0c^2 | dL1/dt INT | dL1/dt SEAM | dL1/dt VTX |
| --- | --- | --- | --- | --- | --- |
| 0.15 | 5.519e-04 | 1.00 | 2.077e-06 | 1.669e-05 | 4.748e-05 |
| 0.30 | 2.337e-03 | 1.06 | 2.215e-05 | 3.548e-05 | 4.652e-05 |
| 0.45 | 5.895e-03 | 1.19 | 3.714e-05 | 3.688e-05 | 3.273e-05 |

`max|v|` scales as **b0c^2** to within 20 % over a 3x change of amplitude — the signature
of the *Maxwell* truncation error, which is what `tests_cs_vertex_cell` section 3 already
isolated (with the field off, the same test gives 1e-17). Halving `b0c` cuts the interior
growth rate **10.7x** and the seam rate 2.1x. So the growing part of `L1(B)` is driven by
the magnetic channel and by nothing else; the hydro half of the update contributes nothing
to it. (The vertex rate is flat in `b0c` because at 0.3 and 0.45 it has already saturated
inside the window.)

## 6. Spherical polar: the same test on a grid with no vertices

`sp_test` iprob 12 is the same toroidal loop (`B = b0c R phihat`, `p0 - b0c^2 R^2`) on a
spherical-polar mesh; run ideal (`eta_ohm_const = 1e-12`), `nx1 x nx2 x nx3 = 16 x 32 x 64`,
t = 2.0. `L1(B)` from its own `.hst`:

| t | 0.063 | 0.126 | 0.251 | 0.500 | 1.001 | 2.000 |
| --- | --- | --- | --- | --- | --- | --- |
| sp `L1-b` | 3.280e-05 | 4.793e-05 | 6.159e-05 | 7.125e-05 | 7.151e-05 | 7.182e-05 |
| sp `L1-v` | 4.850e-05 | 8.979e-05 | 1.449e-04 | 1.820e-04 | 1.741e-04 | 2.472e-04 |

**The sp error SATURATES by t ~ 0.5 and is then flat to t = 2.0** (+0.8 % over 1.5 time
units, against +130 % for the cs interior over the same window). But this is a weaker
control than it looks, and the honest reading is in the finaliser's own split:

```
L1 (per unit volume): B(faces) 1.427e-04     POLAR rows: L1 B 2.595e-03 (1.92 % of volume)
                                             INTERIOR:   L1 B 8.243e-06
```

On spherical polar the toroidal loop is **grid-aligned** — `B_phi` is a coordinate
component and `B_r = B_theta = 0` exactly — so the interior has almost no angular truncation
error to accumulate (8.2e-06, 17x below the cs n = 16 interior), and what error there is
sits entirely in the polar rows. The sp run therefore shows that *this* problem is nearly
exactly representable on sp, not that the cs growth is cs-pathological. **The decisive
evidence remains the h^2 convergence of the cs growth rate in section 1**, which is a
property of the cs scheme measured on the cs grid and needs no cross-grid comparison.

## 7. Verdict

**(i).** The post-t~0.2 growth of `L1(B)` on the cs field loop is the scheme's own
second-order error accumulating secularly, not an instability and not a vertex-driven
drift.

1. The exact solution is **static**, so `L1` against the initial field IS the error; there
   is no advected reference to argue about.
2. The mechanism is a **steady spurious velocity field** of Linf amplitude ~1.17e-3 at
   n = 64 (v ~ h; `KE ~ h^3`), established in the first few steps and then constant to
   0.5 % for the rest of the run, which advects and shears the loop. `ME` is conserved to
   1e-5 and `div B` to 1e-14, so this is displacement, not dissipation of the loop.
3. It is **not exponential**: the interior grows close to linearly in t, the vertex band
   decelerates by 6x over the run, `KE` is flat, and `dln KE/dt` *falls* with refinement
   (0.051 / 0.027 / 0.014).
4. **`dL1/dt` converges at 1.88-2.14 in all three regions** over two refinement steps, and
   the error at the fixed time t = 2.0 is a clean 1.88-2.02. An instability's rate does not
   converge; a 2nd-order truncation drift's does, and this one does.
5. **(iii) is refuted**: `max|v|` is 21-22 % *smaller* at the cube vertex than in the
   interior at both n = 32 and n = 64, the vertex is where the error saturates first, and
   turning the vertex halo fills on or off moves the rate by under 1 %. The vertex's 2.1x
   larger rate is the same O(h^2) Maxwell-stress *constant* already attributed in
   `tests_cs_vertex_cell` section 4, not a separate mechanism.
6. **No knob is needed and none cures it.** `ppmx`/`wenoz` cut the interior rate 1.5x and
   the residual KE drift 4x, and `cs_wellbalanced_src` buys ~12 % uniformly; `fofc`,
   `cs_lowbeta_fallback`, `cs_vertex_fill` and `cs_vertex_fill_cc` do nothing at all.
7. The one **real** instability on this test is behind `<mhd>/bs_emf = true` (the plain
   four-face corner EMF): `dln KE/dt = 1.27`, `KE` x7.8, `L1` 10-20x worse. The shipped
   default (GS05 upwind corner EMF) is what prevents it. That is worth knowing and worth
   leaving alone.

Recommendation: **close gate (c)**. If a future test wants a tighter number on this
problem, `wenoz` is the lever, and `tlim` must be held fixed across resolutions (the
`tlim = 0.25` sampling artefact of `tests_cs_regions` is the same trap).

## Source change

One diagnostic only, in the allowed file:
`src/pgen/cs_test.cpp` — `CSTestLoopHistory` (a `user_hist_func` for iprob 11) plus the
`cs_hist_nband` file-scope copy of `<problem>/conv_nband`. **No solver file was touched**,
and nothing is committed. It is validated against `CSTestResistCheck` to every printed
digit (section "How it was measured").

## Clean-up

Already done: `out/`, the two binaries and both build directories (`build_csloop`,
`build_csloop_sp`, 123 MB) are deleted. What is kept is 144 KB: `README.md`, `growth.py`,
the two inputs, the three run scripts, `logs/` (the full solver logs, including the
`CSTestResistCheck` tables), and `hst/` (the `*.user.hst` time series every table above is
built from, so `growth.py hst/*.user.hst` reproduces them with no rerun).

The sp control at doubled resolution (`sp_n32`, 32 x 64 x 128) was **cancelled** at
t = 0.875 and is not reported: section 6's conclusion rests on the finaliser's
interior/polar split, which n = 16 already gives, and the decisive evidence is section 1.
