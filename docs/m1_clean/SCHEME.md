# The clean M1 radiation scheme: method and properties

This page describes the grey radiation-hydrodynamics scheme selected by
`<rad_m1>/transport = implicit` with `scheme = clean` (the default with implicit transport).
It stays at the level of the method and its properties. The input keys are in
[KEYS.md](KEYS.md), and the tests named below are in `tests/m1_clean/` (see its README for
the decks and the exact criteria).

## 1. The method in brief

**Equations.** Two-moment (M1-type) grey radiation transport for the radiation energy
density $E$ and flux $\mathbf F$, in the mixed frame: lab-frame moments with
comoving-frame opacities, keeping the $O(v/c)$ terms and the $O(\beta^2)$ terms that are
leading order in dynamic diffusion (Krumholz et al. 2007; the β terms follow Ma et al. 2025,
eqs. 56-57):

$$
\partial_t E + \nabla\!\cdot\mathbf F = -c\,G^0,\qquad
\partial_t \mathbf F + c^2\,\nabla\!\cdot(\mathsf D\,E) = -c^2\,\mathbf G ,
$$

where $G^0$ and $\mathbf G$ are the energy and momentum exchange with the gas
(absorption, emission, scattering). The gas receives exactly the opposite:
$\partial_t(\rho e_{\rm tot}) \mathrel{+}= c\,G^0$, $\partial_t(\rho\mathbf v) \mathrel{+}= \mathbf G$.
The gas EOS is gas-only (any AthenaK EOS, including tabulated ones); radiation pressure acts on
the gas only through $\mathbf G$. The scheme requires the full speed of light, $\hat c = c$.

**Closure: variable Eddington tensor from rays.** $\mathsf D = \mathsf K/J$ comes from a
formal solution of the grey transfer equation $\hat n\cdot\nabla I = \chi(S - I)$ along a
fixed set of global directions (the `vet_gd` rays: a Gauss-Legendre product set,
level-symmetric, geodesic or HEALPix set). The sweep is a short-characteristics sweep over
the whole mesh, shell by shell on spherical-polar and cubed-sphere meshes (Cartesian meshes
use the `vet_sc` short-characteristics sweep, optionally with the same sets). The source
function is the state's own ($\epsilon\,aT^4 + (1-\epsilon)E$, or an absorption/scattering
split), so the closure does not depend on $\Delta t$. Besides $\mathsf D$, the rays give the
half-range fluxes (outgoing and incoming parts of the flux) and the rays' own flux through
each face, which the face fluxes use in optically thin regions. Normally a step costs one ray
sweep: the second stage uses an extrapolation of the closure in time.

**Staggered layout.** $E$ lives at cell centres; the radiation flux lives on the cell faces,
one normal component per face. The face flux of a stage is a blend of

- a compact diffusion / wave flux, $-(c/\chi_f)\,\partial_n(\mathsf D E)$ in the optically
  thick limit and the staggered wave operator in the thin limit (no added dissipation needed);
- a kinetic (upwind discrete-ordinates) flux from the rays' half ranges, which takes over where
  the field is anisotropic;
- the advective flux $\mathbf a E$ of radiation carried by the gas, and lagged terms (the
  off-diagonal stress and a consistency correction towards the rays' own thin-limit flux).

**Time integration.** IMEX. The hydrodynamics is AthenaK's explicit Heun (`rk2`) integrator,
unchanged. Radiation transport and the gas-radiation exchange are implicit, with Alexander's
two-stage L-stable SDIRK tableau (`time2_tableau = sdirk2`); each stage sees the gas at the
same time as the explicit part, so the IMEX scheme is second order. Backward Euler is the
fallback (first step, or a step whose stage failed).

**Implicit stage solve.** In each stage the face momentum equation is eliminated face by face,
which leaves a scalar system for $E$. It is solved by Newton passes: each pass is a
right-preconditioned BiCGStab with a line-relaxation multigrid preconditioner, followed by an
exact per-cell root for the gas temperature (any EOS), then a line solve that recovers the face
fluxes. Convergence is measured in each cell's own energy budget of the stage, with the
tolerance tied to the step's time-error estimate.

**Well-balanced formulation.** For a hydrostatic envelope in radiative equilibrium, the
problem generator supplies the reference radiative acceleration of the initial state; the
well-balanced hydro already holds it, and the radiation hands the gas only the residual
force. The face flux equations are written as deviations from the rays' own steady state,
and a steady start (`clean_init = steady`) relaxes the initial radiation to the scheme's own
discrete steady solution in the given gas.

## 2. Properties

<!-- SKELETON: the wording of each property (exact / booked / approximate) and its
mechanism are being audited and will be filled in from PROPERTIES.md. -->

### Energy conservation (gas + radiation)

- Test: `conservation_sp`, `conservation_cs` (dE/E <= 1e-13).
- Mechanism: TBD (PROPERTIES.md).

### Momentum

- Test: none in the suite (the conservation tests check energy, L_z and mass).
- Mechanism: TBD (PROPERTIES.md).

### Angular momentum

- Test: `conservation_sp`, `conservation_cs` (drift of the history column `Lz_tot` <= 1e-13; `Lz_tot` = gas + radiation + the L_z ledger; on the cubed sphere the radiation part is that of the face momenta).
- Mechanism: TBD (PROPERTIES.md).

### Mass

- Test: `conservation_sp`, `conservation_cs` (dM/M <= 1e-13).
- Mechanism: TBD (PROPERTIES.md).

### Positivity by construction

- Test: realisability counts printed every solve; no dedicated test.
- Mechanism: TBD (PROPERTIES.md).

### Well-balanced hydrostatic and radiative equilibrium

- Test: `wb_sp`, `wb_cs` (max |v|/c_s <= 1e-8 over the run).
- Mechanism: TBD (PROPERTIES.md).

### Second order in space and time

- Test: `order_space` (linear radiation-hydro wave, Cartesian; order of the finest pair >= 1.7), `order_time` (deck D, dt halved on a fixed grid; second order is reached as dt -> 0: successive orders 1.57, 1.61, 1.84).
- Mechanism: TBD (PROPERTIES.md).

### Bitwise restart and dt retry

- Test: `restart_sp`, `restart_cs`, `dt_retry` (difference 0; a retry with factor 1 repeats the step exactly).
- Mechanism: TBD (PROPERTIES.md).

### Reproducibility across rank counts and MeshBlock layouts

- Test: `layout_tile` (tile mode: 1 x 2 and 2 x 2 MeshBlocks vs one MeshBlock; on this release they differ by about 1e-10 of the largest value after 8 cycles, fixed in the next release). Rank-count independence (block mode): not in the suite.
- Mechanism: TBD (PROPERTIES.md).

### CPU / GPU consistency

- Test: the suite run on a GPU binary (`tests/m1_clean/README.md`, "GPU builds"); see [BUILD.md](BUILD.md) for a known toolchain problem.
- Mechanism: TBD (PROPERTIES.md).

### No silent fallbacks

- Test: unsupported options stop with a FATAL error naming the key; handled events print `CLEAN WARNING` lines ([OUTPUTS.md](OUTPUTS.md)); a healthy run prints none.
- Mechanism: TBD (PROPERTIES.md).

### Cartesian, spherical-polar and cubed-sphere meshes, one implementation

- Test: the `*_sp` and `*_cs` variants above; `order_space` on a Cartesian mesh. Uniform meshes only (no SMR/AMR).
- Mechanism: TBD (PROPERTIES.md).

## 3. Code map

Paths relative to `src/`.

| part | files |
|---|---|
| key parsing, defaults, FATAL checks | `rad_m1/rad_m1.cpp`, `rad_m1/rad_m1_clean.cpp` (`CleanPreflight`, `CleanInit`) |
| stage solve, face fluxes, write-back | `rad_m1/rad_m1_clean.cpp`, `rad_m1/rad_m1_clean.hpp` (`CleanSolve`, `CleanFaceFlux`, `CleanWriteBack`) |
| linear solver and preconditioner | `rad_m1/rad_m1_clean_solver.cpp` (`CleanBiCGStab`, `CleanStagVCycle`) |
| closure from the rays | `rad_m1/rad_m1_clean_closure.cpp` |
| ray sweeps and direction sets | `rad_m1/rad_m1_vetgd.cpp` (spherical polar), `rad_m1/rad_m1_vetgd_cs.cpp` (cubed sphere), `rad_m1/rad_m1_vet.cpp` (Cartesian) |
| IMEX stages, restart state | `rad_m1/rad_m1_time2.cpp`, `utils/imex_sdirk2.hpp` |
| task lists | `rad_m1/rad_m1_tasks.cpp`, `driver/driver.cpp` |
| dt retry checkpoint | `driver/dt_retry.cpp` |
| layout-invariant sums | `utils/repro_sum.hpp`, `rad_m1/rad_m1_layreduce.hpp` |

## 4. How to cite

- AthenaK: J. M. Stone et al. (2024), arXiv:2409.16053.
- This radiation scheme: a method paper is in preparation; until it appears, cite this
  document with the commit hash of the code you ran.
- Methods it builds on: Ma, Pakmor, Justham & de Mink (2025), arXiv:2503.16627 (mixed-frame
  β terms); Jiang (2021), ApJS 253, 49 (implicit discrete-ordinates radiation, level-symmetric
  sets, scattering split); Davis, Stone & Jiang (2012), ApJS 199, 9 (VET from short
  characteristics); Krumholz et al. (2007), ApJ 667, 626 (mixed-frame sources); Alexander
  (1977), SIAM J. Numer. Anal. 14, 1006 (SDIRK2); Zalesak (1979), J. Comput. Phys. 31, 335
  (flux-corrected transport).
