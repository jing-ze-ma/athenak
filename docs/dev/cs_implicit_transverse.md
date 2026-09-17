# The implicit transverse radiative operator on the cubed sphere

`<hydro>/rad_implicit_ang` (`src/diffusion/conduction_transverse.cpp`, coefficients in
`Conduction::BuildAngularCoeffs`, `src/diffusion/conduction.cpp`) is refused on any
curvilinear mesh today:

```
if (pp->pmesh->use_cubed_sphere || pp->pmesh->use_spherical_polar) {
  ... "rad_implicit_ang is Cartesian-only in this version"
```

This note is the design for lifting the **cubed-sphere** half of that guard.  Spherical
polar keeps its fatal: its pole rows are a separate problem (the `x2` lines end on the
axis, where the operator is singular and the existing polar machinery averages rows).

---

## 1. What the operator is today (Cartesian)

With `e` the internal energy density and `T` linearised about the frozen pre-update
state, `T_i = T*_i + alpha_i (e_i - e*_i)`, `alpha_i = 1/(rho_i c_{v,i})`, the increment
`y_i = e_i - e*_i` obeys

```
dy_i/dt = M(y)_i = (1/V_i) sum_f s_f C_f (T_j - T_i),      y_i(0) = 0,
```

integrated over `t = 0 -> beta_dt`.  `C_f = A_f K_f/dl_f` is the frozen face conductance
built by `BuildAngularCoeffs` from the same `RadFaceKappa` x flux limiter x tau-blend
weight the explicit `x2`/`x3` face flux uses.  On a Cartesian mesh the code exploits
`A_f = 1/dx`, `V_i = 1` (the flux-divergence form the RK update applies), so `C_f` is
`K_f/dx^2` and **no volume appears anywhere in `conduction_transverse.cpp`**.

Two solvers integrate the same `M`:

* `rad_ang_solver = sts` (default): an RKL1 super-time-stepping loop, `s` substages of
  the 5-point stencil, `s` from a global Gershgorin radius.  It needs only a halo
  exchange of the increment; there is no line structure at all.
* `rad_ang_solver = adi`: one alternating-direction step.  `M(y) = b + A y`,
  `A = A2 + A3`; a sequence of backward-Euler sub-steps
  `(I - dts A2) Y1 = y_in + dts b2`, `(I - dts A3) y = Y1 + dts b3`
  (`rad_adi_scheme = lod`), Richardson-extrapolated in the step size for second order
  (`lod2`, the production scheme).  Each sweep is a **tridiagonal line solve that crosses
  MeshBlocks and ranks**, solved by partition: one Thomas factorisation with three
  right-hand sides per block piece, a 2-unknowns-per-block reduced interface system
  gathered around the ring of blocks of the line and solved redundantly, then local back
  substitution.  The ring is cyclic, which is why a direction split over more than one
  MeshBlock is required to be periodic.

Both are written in **flux form**: the two cells sharing a face form the identical
product `C_f (Th_j - Th_i)` from identical operands, so `sum_i V_i y_i = 0` to round-off.

A face on a physical (non-`block`, non-`periodic`) `x2`/`x3` boundary is **closed**:
`open2`/`open3` test `mb_bcs` and return false, and the face carries nothing.  There is
no Dirichlet condition -- the increment is simply not defined outside the mesh.

---

## 2. The cubed-sphere operator

`x1` is the radius; `x2 = xi` and `x3 = eta` are the gnomonic panel coordinates in
`[-1,1]`, six panels.  At fixed `r` the two coordinate lines meet at an angle `alpha`
with (`Coordinates::sin_cell`, `Coordinates::cos_cell`, both `(m,k,j)`)

```
sin alpha = sqrt(1 + X^2 + Y^2)/(C D),    cos alpha = -X Y/(C D),
X = tan(xi), Y = tan(eta),  C = sqrt(1+X^2), D = sqrt(1+Y^2),
```

so the surface metric in the arc-length frame is `g = [[1, cos a],[cos a, 1]]`, with
inverse `g^{-1} = [[1, -cos a],[-cos a, 1]]/sin^2 a`.  **`|cos alpha| <= 1/2`**, with the
maximum at a cube vertex (`X = Y = 1`): the cross term is always strictly bounded by the
diagonal ones, `|g^{23}| < sqrt(g^{22} g^{33})`.

The face-normal derivative on a `xi`-face is therefore not `dT/dl_xi`:

```
dT/dn|_xi  = ( dT/dl_xi  - cos a  dT/dl_eta ) / sin a
dT/dn|_eta = ( dT/dl_eta - cos a  dT/dl_xi  ) / sin a
```

which is exactly what the **explicit** operator already does
(`AddIsotropicHeatFluxRadiative`, kernels `radcond2`/`radcond3`, gated on
`use_cubed_sphere && rad_cs_exact`) and what the frozen coefficients already carry in
part.  Writing the total flux through an `x2` face (area `A_f = area.x2f`, arc length
`dl_f = (dx2(j-1) + dx2(j))/2`, `K_f` the limited face conductivity in code units):

```
Phi_f = A_f * ( -K_f * dT/dn )
      = -[A_f K_f/(dl_f sin a)] (T_j - T_{j-1})
        + [A_f K_f/(dl_f sin a)] * dl_f * cos a * (dT/dl_eta)_f
      = -C_f (T_j - T_{j-1}) + G_f * ge_f,
```

with

```
C_f = A_f K_f/(dl_f sin a)      <-- ALREADY what cap_c2/cap_c3 hold (curvilinear branch)
G_f = C_f * dl_f * cos a        <-- NEW
ge_f = the eta-derivative at the face, the two-cell average the explicit kernel uses:
       0.5*[ (T(k+1,j-1)-T(k-1,j-1))/Dz(k,j-1) + (T(k+1,j)-T(k-1,j))/Dz(k,j) ],
       Dz(k,j) = 0.5 dx3(k-1,j) + dx3(k,j) + 0.5 dx3(k+1,j).
```

So **`BuildAngularCoeffs` already builds the right `C_f` on the cubed sphere**
(`af = area2_`, `/(dl*sn)`), and already evaluates the cross term -- but only inside the
flux limiter's `gradn` argument.  What is missing from the implicit operator is:

1. `G_f ge_f`, the cross flux itself;
2. the cell volume: `M(y)_i = (1/V_i) sum_f s_f Phi_f` with `V_i = volume(m,k,j,i)`,
   not `1`;
3. the `panel` boundary flag: `open2`/`open3` close a panel seam today, so a cubed-sphere
   run would silently lose every cross-panel transverse flux;
4. the same three things in the Gershgorin radius and in the conservation diagnostic
   (which divides by `dx1 dx2 dx3` today).

### The split: diagonal implicit, cross explicit

`ge_f` couples `j` to `k`, i.e. the stencil is 9-point (the four diagonal neighbours
enter).  A tridiagonal sweep cannot carry it, and neither can the RKL1 M-matrix argument.
The standard treatment, and the one adopted here:

* the **diagonal** parts `C_f (T_j - T_i)` stay exactly where they are -- in the implicit
  sweeps of the ADI, in the stencil of RKL1;
* the **cross** parts `G_f ge_f` are evaluated **explicitly**, from a lagged temperature,
  and enter as a source.

For RKL1 the lag is zero: `M(y)` is evaluated at the current substage value
`Th = T* + alpha y`, cross term included, so RKL1 integrates the *exact* cubed-sphere
operator; only the Gershgorin bound has to be widened to cover the extra entries.

For ADI the cross term is part of the affine source `b = M(0)`, i.e. it is frozen at
`T*` over the stage.  `rad_adi_cross_iter > 1` re-evaluates it at `Th = T* + alpha y` and
repeats the sweeps, a fixed-point iteration whose contraction factor is the ratio below.

### Why the explicit cross term is stable

For one Fourier mode on the panel, with `A_d` the diagonal part and `A_c` the cross part
and `alpha = V = 1`, write `A_d T = -a T`, `A_c T = -c T`.  The LOD/backward-Euler step
with the cross term in the source gives

```
T_new/T* = 1 - tau (a + c)/(1 + tau a) = (1 - tau c)/(1 + tau a).
```

Stability needs `|1 - tau c| <= 1 + tau a`, i.e. `-2/tau <= c <= a + 2/tau`, and the
symbols supply `|c| <= a/2`:

```
a = [ (2 - 2 cos th2) + (2 - 2 cos th3) ] / h^2      (diagonal, g^22 = g^33 = 1)
c = 2 g^23 sin th2 sin th3 / h^2,    |2 g^23| <= 1
|c|/a <= sin th2 sin th3 / [ (2-2cos th2) + (2-2cos th3) ] <= (1 + cos th)/4 <= 1/2.
```

Two things follow, and they are the honest characterisation of the scheme:

* it is **unconditionally stable** -- `|T_new/T*| <= max(1, (tau a/2 - 1)/(1 + tau a)) < 1`
  for every mode and every `tau`;
* the *worst* ratio `|c|/a = 1/2` is reached as `th -> 0`, i.e. at the **longest**
  wavelengths, where `tau a` is small and the explicit treatment is accurate anyway; at
  the Nyquist checkerboard, which is what the implicit operator exists to damp, the
  central cross difference has symbol `sin(pi) = 0` and the cross term **vanishes**.  The
  residual worst case is the mid-band `th2 = th3 = pi/2`: `|c|/a = 1/4`, amplification
  `-1/4` per stage instead of the exact `~0` -- damped by 4x per stage, but with a sign
  flip, so a mid-band mode rings down rather than being annihilated.  `rad_adi_cross_iter`
  is the switch that removes it (contraction 1/4 per iteration).

### Conservation

Take the operator in the form

```
V_i dy_i/dt = sum_f s_f [ -C_f (Th_j - Th_i) + G_f ge_f ] * (-1)^{...}
```

i.e. every face contributes one number `Phi_f` to the cell on each of its sides, with
opposite signs.  `C_f`, `G_f`, `dl_f`, `cos a` and the stencil of `ge_f` are all built
from operands both cells can see (cell-centred arrays plus their ghost ring), and both
cells evaluate the *same expression*, so `Phi_f` is bitwise equal and opposite and

```
sum_i V_i y_i = (flux through the physical x1/x2/x3 boundaries alone)
```

to round-off.  This is exactly the property the Cartesian code has, carried over by
writing the volume division *outside* the face sum and never inside it.  The one place
the equality is not bitwise is a **panel seam** (section 4).

---

## 3. The ADI line structure on the cubed sphere

`x2` lines do not continue across a panel seam: the neighbouring panel's `x2` may be the
first panel's `x3` (signed axis swap, `PanelEdgeMap`), and the two charts' cells do not
even coincide along the seam -- the halo is filled by a quadratic **along-seam resample**
(`bvals_cc.cpp`, `cs_seam`).  There is no cyclic ring of blocks to gather.

Therefore, on the cubed sphere:

```
linked face  (in the tridiagonal line)   : BoundaryFlag::block, periodic
open face    (carries flux)              : block, periodic, AND panel
closed face  (carries nothing)           : everything else (reflect, outflow, user, ...)
```

The reduced interface system is then assembled on a **chain** instead of a ring: the
existing cyclic assembly already degenerates correctly when the end couplings `a_1`,
`c_n` are zero, which is what a non-linked face gives it, so the ring gather and the
reduced solve need no change beyond the block count per line becoming *per panel*
(`mesh_indcs.nx2/mb_indcs.nx2` within one panel -- the panel trees are separate, so this
is already the per-panel count) and the periodicity check being replaced by "periodic or
panel".

## 4. The seam sub-step

A panel face left as a pure source term would be an explicit face embedded in an implicit
operator: for the mode that is uniform on each side and jumps at the seam, every other
face carries no flux, the amplification is `1 - 2 z` with `z ~ 10^3` in production, and
the operator blows up.  (The RKL1 path has no such problem: it is stable for every face
in its stencil by construction, and a panel face is just another open face to it.)

The ADI path therefore closes panel faces during its sweeps and applies them in a
separate, conservative, unconditionally stable **pair-implicit sub-step** afterwards.
For each panel face `f` with owned cell `A` and halo cell `B`:

```
Th_A = T*_A + alpha_A y_A,    Th_B = T*_B + alpha_B y_B          (y halo exchanged)
delta = tau [ C_f (Th_B - Th_A) + s_f G_f ge_f ]
        / ( 1 + tau C_f (alpha_A/V_A + alpha_B/V_B) )
y_A += delta/V_A
```

This is the exact backward-Euler solve of the isolated two-cell system.  It is

* **L-stable**: `delta -> (Th_B - Th_A)/(alpha_A/V_A + alpha_B/V_B)` as `tau -> inf`,
  i.e. the pair is relaxed to its common temperature and no further -- no overshoot for
  any `tau`;
* **conservative**: `B` forms `delta_B` from the same operands with the roles swapped and
  gets `-delta_A`, up to the along-seam resample.

The residual seam non-conservation is the resample error, and it is **the same error the
explicit operator already makes**: `radcond2`/`radcond3` also form each seam flux twice,
once per panel, from resampled halo temperatures.  The test plan measures it against that
reference rather than against zero.

## 5. Halo requirements

* `ge_f`/`gx_f` read the **diagonal** `(j+-1, k+-1)` ghosts, so
  `rad_tr_halo_faces_only` (which drops the `x2x3` edge/corner buffers) must be refused on
  the cubed sphere.
* `T*` and `alpha` are already exchanged through `pbval_tr`
  (`MeshBoundaryValuesCC`), which handles panel seams for a cell-centred scalar: the
  momentum transform is gated on `v == IVY || v == IVZ` and `pbval_tr` carries one
  variable, so only the index permutation and the along-seam resample apply.
* The ADI seam sub-step needs one extra exchange of `y` per sub-step.
* A **cube vertex** corner ghost is not exchanged at all (`IsCubeVertexCorner`); it is
  filled by `mesh/cs_vertex_fill`.  The cross-term stencil reads it, exactly as the
  explicit `radcond2` does -- same exposure, no new one.

## 6. What stays refused

* spherical polar (`use_spherical_polar`) -- unchanged fatal;
* SMR/AMR -- unchanged fatal (the increment exchange would have to prolongate/restrict);
* `rad_sts_all` on a curvilinear mesh -- unchanged (one global substage count would pay
  the radial stiffness on every transverse face);
* `rad_tr_halo_faces_only` on the cubed sphere -- new fatal.

## 7. Bitwise requirement

Every term added here is guarded by `use_cubed_sphere` (and, for the cross term, by
`rad_cs_exact`, exactly as the explicit operator is).  The Cartesian arithmetic -- the
order of the products in the face flux, the absence of a volume division, the `open2`
truth table for `block`/`periodic` -- is not touched, so a Cartesian run is bitwise what
it was.  The regression test (a) below is the proof obligation.

---

## 8. Test plan

**(a) Cartesian bitwise regression.**  `inputs/tests/rad_transverse_gauss.athinput`
(pgen `rad_diff2d`) and the He-star box smoke `he_box_w8.athinput` at `nlim = 50`, run
with the head binary and the new one for `rad_ang_solver = sts` and `= adi`
(`rad_adi_scheme = lod2`).  `cmp` of the binary dumps must be identical.

**(b) Cubed-sphere diffusion.**  `cs_test` `iprob = 15` already sets a uniform density at
rest with `T = T0 (1 + amp f)` and `f` an `l = 2` combination of direction cosines, no
radial dependence, so the radial fluxes vanish identically and the transverse operator is
the only one that moves anything; `inputs/tests/cubed_sphere_raddiff.athinput` is the
input.  Extended with an `l = 6` zonal harmonic `P_6(cos)` (`lharm = 6`), for which
`Laplacian_sphere f = -42 f/r^2`.  A harmonic of degree `l` decays as

```
A(t)/A(0) = exp( -l(l+1) D t / r^2 ),    D = K/(rho c_v),
```

measured by the volume-weighted projection `A(t) = sum_i V_i (T_i - T0) f_i / sum_i V_i f_i^2`.
Reported for
  * the new implicit operator, `rad_ang_solver = sts` and `= adi` (`lod2`), and
  * the existing **explicit** cubed-sphere operator (`rad_implicit_ang = false`, i.e. the
    reference `radcond2`/`radcond3` path, with `rad_cap_ang = 0`),
at 32 and 64 cells per panel edge.  Targets: rate error < 2 % at 32, second order between
the two resolutions; shell-integrated energy `sum_i V_i (e_i - e_i(0))` conserved to the
level the explicit reference achieves; the seam signature reported as
`max |T - T_fit|` over the cells adjacent to a panel edge against the same over the
panel interiors.

**(c) Stability.**  The same test at a transverse diffusion number `D dt/dx^2 ~ 10^3`
(raise `rad_kappa_fac`, or run at fixed `dt` with the hydro CFL doing the limiting): the
solution must stay bounded and monotone (no sign change of `A(t)`, no growth).
