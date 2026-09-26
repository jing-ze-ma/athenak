# Task for viper: implicit M1/VET on a REAL spherical wedge (rad-hydro) — handed over from Caltech 09-25 PDT

**User 09-25 (Caltech):** "try to have implicit vet working on a real wedge with rad hydro, and do the same
as you have done on the box convection". This settles open decision #2 of HANDOVER-2026-09-26 (section 6):
build a gravity-bearing spherical M1 wedge. It is handed to **viper** because Caltech is GPU-bound (the H200
debug queue allows one job at a time). No code was written on Caltech; the design research is summarised
below.

Scope: a **test-size** development and performance bed, not a production model. Work on a new branch
(e.g. `m1-wedge`) from `rt-integration` and merge when gated, as usual.

## What exists today

- The only sp M1 setup is test **T-S4**: `rad_m1_beam`, `m1_test = sph_atm`. It is a grey atmosphere with
  **no gravity and the gas held fixed**. All the fast5-sp / vet_col / cfl timings used it.
- The 4 Msun presupernova He-star plan (memory `he4-presn-global-plan`) used **grey two-stream RT** on the
  cubed sphere, not M1. The He star parameters (M 3.15, log L 4.78, Teff 49 kK, R 2.37e11 cm; FeCZ ~0.64–0.97
  R, tau 11–208) and the spherical 1-D structure work (`column_sph.py`, the ic_profile route) can be reused.
- The Cartesian He box, `src/pgen/box_convection.cpp`, is the physics template. Its input is
  `docs/handover/caltech-2026-09-26/inputs/he_box/box.athinput`.

## box_convection facts that matter for a sphere (researched on Caltech, file:line from rt-integration 9c4b57f7)

- **It fatals on sp/cs** (`box_convection.cpp:1072-1076`). **x1 is vertical everywhere:**
  - the potential `g0 (z - zmin)`;
  - the IC column;
  - the BCs;
  - the history planes;
  - `wall_noflux`;
  - the horizontal area (x2/x3 are periodic).
- **Gravity** is a constant-g0 user source (`BoxConvSrcs`, :2626-2677). Under well-balancing it is the WB
  background pressure difference `bdt (pr - pl)/dz`, commented "Cartesian: equal face areas" (:2651). The
  sphere needs a point mass `GM/r^2` and area-weighted WB faces.
- **EOS:** general tabulated, gas only. rad_m1 fatals if `eos_radiation = true` (`rad_m1.cpp:463-465`).
- **Opacities are handed over by the pgen:**
  - `ReadOpacityTable` (:775-825) reads Rosseland + Planck on a shared grid.
  - They go to `pm1->SetOpacityTables` (:1635).
  - rad_m1 fatals without them.
  - A sphere pgen must do the same.
- **IC:**
  - `ic_profile` is `z rho eint` in cgs, interpolated in log (:1239-1290).
  - Radiation IC: `m1_ic_file` `z E F` (:2410-2503).
  - WB reference acceleration: `wb_arad_file` `z a_rad` (:1390-1459) and `SetForceReference` (:1695-1745).
  - The data files are `HE_BOX_DATA/`: ic_m1_V3edd_pgen.txt, m1_rad_ic_V3edd.txt, rosseland/planck he
    tables, arad_V3edd.txt.
- **BCs:**
  - Hydro: bc_mode 3, a WB hydrostatic walk at both x1 walls.
  - M1 ghosts: flux bottom (E from the diffusion extrapolation) and a dark top ghost.
  - Implicit solver faces: `implicit_bc_x1min = flux`, `x1max = marshak` (`rad_m1_implicit.cpp:1011-1016`).
  - The pgen moves `rad_flux_inner` into rad_m1 so that the luminosity is not injected twice (:1658-1668).
- **History:** F1 top/mid/bot vs F_bot, V1max, Etot, Fres, KEcol, and `rt_profile.bin`.

## Suggested design (confirm or change)

- A new pgen, or an sp mode of an existing one: a **spherical-polar wedge, point mass, open/WB radial walls**.
  Theta/phi walls are reflecting or periodic.
- Physics as in the box: the same tables, hesdirk2, implicit bicgstab, **vet_col** (sp), and precond `mg_gc`
  with levels 1 (the sp default).
- **IC** from a 1-D spherical hydrostatic + radiative profile (the `ic_profile` route; reuse the he4
  spherical column tools).
- **Region:** below the photosphere down to tau ~100–300 around the FeCZ.
- **Test grid:** about 64–128 radial x 32 x 32 cells, a 10–20 degree wedge.
- **Time scheme:** reference cfl 0.3; cfl 0.9 as an extra arm (the He-box KE question is open).
- **Accuracy region** (user rule): judge only below the photosphere, tau >~ 1, excluding the top cells and the
  sponge. The deviation bar is the round-off spread: 1 vs 2 GPUs, or a 1e-14 kick.

## Gates (same as the box work)

1. CPU:
   - the IC holds or relaxes sensibly (velocities);
   - top luminosity vs the base flux;
   - no NaN/FATAL; NOT-CONVERGED counts;
   - multi-rank;
   - restart bitwise.
2. Existing problems stay **bitwise**: the He box, the dhj WASP-121b gate and the T-S4 wedge, old vs new.
3. GPU (MI300A):
   - 1 vs 2 GPUs;
   - restart bitwise;
   - timing of the sp M1 levers: mg_gc vs rbgs_fwd, vet_col_chunk, op_team_red, predictor, one_pass, EW,
     vet_col_every (relaxation only), cfl 0.9 vs 0.3;
   - a profile.

## Caltech status of the same code (for later GPU comparison on H200)

- The dhj WASP-121b port is merged (`HANDOVER-2026-09-26.md`, top block). The **He box does not yet build under
  nvcc**: `box_convection.cpp:3278` has a local lambda `state_i` captured in the `fill` KOKKOS_LAMBDA.
- A Caltech agent is porting `box_convection` + `rad_m1` on branch `m1-port` (not yet pushed).
- **Any new wedge code should avoid the patterns that break CUDA:**
  - generic or nested helper lambdas inside kernels;
  - host code indexing device Views;
  - DualViews filled on the host without `sync_device()`;
  - class members read inside kernels (capturing `this`).
- Once the wedge is gated on viper, Caltech can build it for H200 and compare.
