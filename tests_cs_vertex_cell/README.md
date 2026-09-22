# Is the CUBE-VERTEX CELL'S OWN UPDATE inconsistent?  Diagnosis, and the answer: no.

`tests_cs_regions` measured a cube-vertex error hot spot (`strat` 2.5-2.9x the interior;
`loop` "order 1.57-1.65 against 2.0" with a VERTEX/INTERIOR ratio growing 1.75 -> 3.06)
and, after excluding the ghost fills, the radial reconstruction, the well-balanced source
and the corner EMF, attributed it to "the incomplete discrete cancellation between the
flux divergence and the geometric source in the skewed cells" — 1 part in 260 at a vertex
against 1 part in 600 in the interior.

This directory takes that last step apart. **Every term is consistent.** The geometry is
exact, the discrete Gauss identity is no worse at a vertex than in the interior, and the
flux-divergence / geometric-source cancellation is EXACT TO ROUND-OFF at a vertex for a
uniform state *and* for the hydrostatic stratification. No fix is shipped, because
nothing was found to fix.

Everything below is serial CPU, Release, `-D PROBLEM=cs_test`, one unmodified binary
(`athena_ref`, built from `fdc08a21` with no source change).

```bash
cmake -B build_csvcell -D PROBLEM=cs_test -D CMAKE_BUILD_TYPE=Release && make -C build_csvcell -j8
cp build_csvcell/src/athena tests_cs_vertex_cell/athena_ref
tests_cs_vertex_cell/run1.sh          # one-step (LTE) per-region scan
tests_cs_vertex_cell/runt.sh          # fixed-time per-region scan (TL, NS, TS in the env)
python3 tests_cs_vertex_cell/gauss.py # the discrete Gauss identity, offline
```

---

## 1. The geometry is in closed form, not quadrature

`src/coordinates/coordinates.cpp` `cscoord` kernel, lines 657-706. Nothing here is a
midpoint rule:

| quantity | line | form |
| --- | --- | --- |
| radial face area | 661 | `r^2 * GnomonicSolidAngle(xl,xr,yl,yr)` — the EXACT solid angle |
| cell volume | 676 | `(1/3)(r_r^3-r_l^3) * ` the same exact solid angle |
| xi / eta face area | 666, 671 | `0.5*(r_r^2-r_l^2) * dth` with `dth = acos(...)` the EXACT arc subtended in that plane — the xi = const surface is a plane through the origin, so a radial sector area is exact |
| edge lengths | 693-706 | `r * ` the SAME arcs the adjoining face areas use |

So there is no quadrature error in areas, volume or edges to find. What *is* a midpoint
choice is the cell centre in the two angles (`CellCenterX` on xi and eta), whereas the
radial centre is the volume centroid (`RadialCentroid`, line 606) — see section 5.

## 2. The discrete Gauss identity: fine, and BETTER at the vertex

`gauss.py` forms `|sum_f sigma_f A_f nhat_f|` exactly as the code does — the code's areas
above, and the face-CENTRE normals the well-balanced source uses — and normalises by the
cell's total surface area. A non-zero value is a cell that "leaks" momentum.

| n | INTERIOR | SEAM | VERTEX |
| --- | --- | --- | --- |
| 16 | 3.179e-05 | 2.956e-05 | **2.285e-05** |
| 32 | 3.960e-06 | 3.698e-06 | **2.921e-06** |
| 64 | 4.940e-07 | 4.625e-07 | **3.698e-07** |

`|sum A n|` itself is O(h^5) (ratios 32.7, 32.4 = 2^5), i.e. O(h^2) once divided by the
volume: a clean second-order consistency error, and at the vertex it is **30 % SMALLER**
than in the interior, at every resolution. The vertex cell does not leak.

## 3. For a uniform state the two halves cancel EXACTLY — also at the vertex

Analytically first. For `v = 0`, `B = 0`, uniform `p` the x2 Riemann flux is
`(F_IM2, F_IM3) = (p, 0)` in the face frame; `GnomonicEquiangleFluxX2`
(`gnomonic_kernels.hpp` 172-186) turns that into

    flx.x2f(IM2) = p * sin_face_xi ,   flx.x2f(IM3) = 0

so the x2 divergence of the xi momentum is `-(A_r s_r - A_l s_l) p / V`. And the default
geometric source (`SrcTermsGnomonicEquiangleImpl`, line 1560) adds `x_ov_rD * p` with

    x_ov_rD = (area2r*sin_face_xir - area2l*sin_face_xil) / volume    (line 760)

— the same number. The cancellation is **structural, not asymptotic**, and the same holds
for `y_ov_rC` on the eta faces and `z_ov_rE` on the radial ones. The *default* cs source
is therefore already well balanced for uniform states; `cs_wellbalanced_src` extends that
to anisotropic states, which is why it changes nothing here.

Measured, on `cs_regions_strat` (iprob 13) at `nlim = 1`, `plm`, varying what is switched
on. `L1(v_t)`, the spurious TANGENTIAL velocity, which the exact solution has none of:

| control | what is in it | n=16 INT / SEAM / **VTX** | n=64 INT / SEAM / **VTX** |
| --- | --- | --- | --- |
| **A** `grav=0 b0c=0` | uniform gas, no field | 1.02e-17 / 8.38e-18 / **2.15e-17** | 1.16e-17 / 1.15e-17 / **1.30e-17** |
| **C** `grav=1 b0c=0` | hydrostatic atmosphere, no field | 8.76e-18 / 8.69e-18 / **1.97e-17** | 1.12e-17 / 1.15e-17 / **1.13e-17** |
| **B** `grav=0 b0c=0.3162` | uniform gas + uniform Cartesian B | 4.37e-06 / 8.47e-06 / **1.42e-05** | 7.72e-08 / 1.66e-07 / **2.89e-07** |
| **D** `grav=1 b0c=0.3162` | the gate's `strat` | 1.61e-05 / 3.15e-05 / **5.30e-05** | 2.73e-07 / 5.89e-07 / **1.03e-06** |

**A and C are round-off, at the cube vertex included.** The hydro half of the vertex cell's
update — flux divergence, geometric source, gravity, face areas, volume — balances
*exactly* there. **100 % of the spurious tangential velocity is the MAGNETIC channel**, and
D is just B with the stratification amplifying it 3.5x (`v = m/rho` and rho falls by e^-5).

## 4. The magnetic channel is SECOND ORDER at the vertex too — a constant, not an order

Control B, one step (so the measured order is spatial order + 1, since `dt ~ h`):

| region | n=16 | n=32 | order | n=64 | order | VTX/INT |
| --- | --- | --- | --- | --- | --- | --- |
| INTERIOR | 4.366e-06 | 6.040e-07 | 2.85 | 7.718e-08 | 2.97 | 1.00 |
| SEAM | 8.472e-06 | 1.223e-06 | 2.79 | 1.655e-07 | 2.89 | 2.15 |
| **VERTEX** | 1.421e-05 | 2.115e-06 | **2.75** | 2.889e-07 | **2.87** | **3.74** |

and the full `strat` (control D) is the same picture: 2.89/2.99 interior against 2.79/2.90
at the vertex. Spatial order **1.87 at the vertex against 1.97 in the interior**. The
vertex carries a 3.3-3.7x larger CONSTANT, and the constant is set by
`|d(tangent basis)/d(angle)|`, which is maximal where the three panels' coordinate lines
meet at 120 degrees. That is the whole of the "1 part in 260 vs 1 part in 600".

Who owns the constant, attributed by switching the two levers that could remove it
(control B, `L1(v_t)` at the VERTEX):

| recon | wb | n=16 | n=32 | n=64 |
| --- | --- | --- | --- | --- |
| dc | off | 1.011e-05 | 1.886e-06 | 3.801e-07 |
| **dc** | **on** | 8.778e-06 | 1.648e-06 | 3.467e-07 |
| plm | off | 1.421e-05 | 2.115e-06 | 2.889e-07 |
| plm | on | 1.499e-05 | 2.128e-06 | 2.860e-07 |
| wenoz | off | 3.911e-06 | 6.166e-07 | 8.714e-08 |
| wenoz | on | 3.259e-06 | 5.072e-07 | 7.200e-08 |

The `dc` + `cs_wellbalanced_src` row is the interesting one. That combination makes the
face states *identical* to the cell state, which is the condition under which the
well-balanced source is documented to cancel the divergence to round-off — and it does not
(8.8e-06, second order). The reason is structural and is **not a defect**: in every sweep
the NORMAL field component is taken from `b0.x*f`, not from the reconstructed `bcc`
(`gnomonic_kernels.hpp` 107-113), so at the two xi faces the flux sees `bn_l` and `bn_r`
while the source sees their cell average — and that difference *is* the discrete
`-(B.grad)B`, a real term of the Lorentz force, not an imbalance. There is no state for
which the Maxwell part can be made to cancel by construction, and its truncation error is
the ordinary O(h^2) one.

## 5. What is left, and what it is not

* Not the ghost fills (already excluded), not the radial reconstruction (already
  excluded), not the areas/volume/edges (section 1), not the Gauss identity (section 2),
  not the hydro cancellation (section 3, round-off).
* The one remaining *choice* that is a midpoint rather than a centroid is the cell centre
  in xi and eta, and the face triads at `(xi_face, eta_centre)` / `(xi_centre, eta_face)`
  (`BuildWBGeometry`, coordinates.cpp 1290-1343) — and the cell-centred field's two
  angular slots, `0.5*(b0.x2f(j) + b0.x2f(j+1))` (coordinates.cpp 1072-1073), where the
  radial slot uses the weighted `CellCenteredRadialFld`. All three are second-order
  accurate at the midpoint and none of them can change the ORDER; at most they move the
  constant, and any of them costs a second-derivative stencil to improve. Not attempted.

## 6. The `loop` "order 1.57" is not a truncation order

This is a retraction of the sharpest claim in `tests_cs_regions/README.md`.

**(a) Most of the reported error is already there at t = 0.** `nlim = 0` — the vector
potential seeding alone, no evolution at all — per-region `L1` of the face field:

| region | n=16 | n=32 | order | n=64 | order | VTX/INT |
| --- | --- | --- | --- | --- | --- | --- |
| INTERIOR | 4.911e-05 | 1.2198e-05 | 2.01 | 3.0384e-06 | 2.01 | 1.00 |
| SEAM | 7.022e-05 | 1.7584e-05 | 2.00 | 4.3983e-06 | 2.00 | 1.43 |
| **VERTEX** | 9.894e-05 | 2.5441e-05 | **1.96** | 6.4362e-06 | **1.98** | 2.01 / 2.09 / 2.12 |

Clean second order in all three regions, with a FIXED vertex ratio. Against the gate's
evolved numbers at `tlim = 0.25` (INT 5.611e-05 / 1.3747e-05 / 3.4410e-06; VTX 9.803e-05 /
3.3124e-05 / 1.0542e-05) the t=0 error is 87 % of the interior error at every n, and
**101 %, 77 %, 61 % of the vertex error** — a share that falls with n, which alone bends
the apparent vertex order downwards.

**(b) The evolved part is a transient, and `tlim = 0.25` samples different phases of it at
different n.** `loop` at n = 32, vertex `L1` against time:

| t | 0 | 0.0625 | 0.125 | 0.25 | 0.5 | 1.0 |
| --- | --- | --- | --- | --- | --- | --- |
| INTERIOR | 1.220e-05 | 1.636e-05 | 1.478e-05 | 1.375e-05 | 1.343e-05 | 2.102e-05 |
| **VERTEX** | 2.544e-05 | 2.677e-05 | 2.538e-05 | 3.312e-05 | 6.389e-05 | 1.002e-04 |

The vertex error *falls* to t = 0.125 and then grows 4x by t = 1.0. At the gate's fixed
`tlim = 0.25`, n = 16 is still on the flat part and n = 64 is already growing, so the
three points that the "order 1.57" is fitted through are not the same solution phase.

**(c) At a time where every resolution is in the same phase the order comes back.** `loop`
at `tlim = 1.0`:

| region | n=16 | n=32 | order | n=64 | order | VTX/INT |
| --- | --- | --- | --- | --- | --- | --- |
| INTERIOR | 8.740e-05 | 2.102e-05 | 2.06 | 5.3027e-06 | 1.99 | 1.00 |
| SEAM | 1.676e-04 | 4.305e-05 | 1.96 | 1.1230e-05 | 1.94 | 1.92 / 2.05 / 2.12 |
| **VERTEX** | 3.551e-04 | 1.002e-04 | **1.83** | 2.5824e-05 | **1.96** | 4.06 / 4.77 / 4.87 |

**1.83 then 1.96, not 1.57 then 1.65** — and the VERTEX/INTERIOR ratio 4.06 -> 4.77 ->
4.87 is *saturating*, not growing. On the same test, at a time where every resolution is
past the transient's onset, the cube vertex is a clean second-order region carrying a
~4.9x constant.

## 6b. Blocking control: 1x1 vs 2x2 MeshBlocks per panel

Control B at n = 32, `meshblock/nx2=nx3=16` (2x2 blocks per panel) against `=32` (1x1):
`L1(v_t)` = 6.0398e-07 / 1.2227e-06 / 2.1150e-06 (INT / SEAM / VTX) in **both**, to every
printed digit. The vertex error does not depend on how the panel is decomposed, which is
one more way of saying it is not a boundary-exchange effect.

## 7. Conclusion

1. **There is no leading-order inconsistency in the cube-vertex cell's update.** Areas,
   volume and edge lengths are exact closed forms; the discrete Gauss identity is O(h^2)
   and 30 % *better* at the vertex; the flux-divergence / geometric-source cancellation is
   exact to round-off at the vertex for a uniform state and for the hydrostatic
   atmosphere. Nothing is shipped, and no switch was added.
2. **The 1:260 is the O(h^2) truncation error of the MAXWELL stress**, and only of the
   Maxwell stress: with the field off, the same test returns 1e-17 at the vertex. Its
   constant is 3.3-3.7x the interior because the tangent basis rotates fastest where the
   panels meet at 120 degrees. It is second order (spatial 1.87 vertex vs 1.97 interior),
   so it is a constant, not an order defect — which is exactly why only the
   reconstruction (`ppmx`/`wenoz`, 3-4x) moves it and the source rebuild does not.
3. **The `loop` vertex order of 1.57-1.65 is an artefact** of (i) a t = 0 seeding error
   that is 61-101 % of the total and second order everywhere, and (ii) a transient whose
   onset time depends on n, sampled at a fixed `tlim`. At `tlim = 1.0` the vertex order is
   1.83 (n 16->32) and 1.96 (n 32->64), and the vertex/interior ratio SATURATES at 4.9
   instead of growing. The vertex is a 2-5x error hot spot on this test; it is not a
   lower-order region.
4. Open, and not addressed here: the slow growth of the `loop` error after t ~ 0.2 (all
   three regions grow together at late times, vertex-heaviest). That is a stability
   question, not a consistency one, and it is what the gate's fixed-`tlim` order was
   actually measuring.

## Clean-up

```bash
rm -rf tests_cs_vertex_cell/out tests_cs_vertex_cell/athena_ref build_csvcell
```
Keep `logs/`, `README.md`, `gauss.py`, `run1.sh`, `runt.sh`.
