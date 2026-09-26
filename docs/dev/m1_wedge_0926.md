# m1-wedge (09-26): implicit M1 + vet_col on a gravity-bearing spherical-polar rad-hydro wedge

Task: `docs/handover/TASK-2026-09-26-m1-real-wedge.md` (rt-integration 206318c9). Branch `m1-wedge`
from `c1dc8226` (the base of every bitwise gate below). Run tree, scripts, inputs and logs:
`/viper/ptmp2/jinma/sprhd_0926` (RESUME.md there).

## 1. Audit: what the sp M1 path supports when the gas moves (file:line at c1dc8226)

Verified claims of the TASK appendix, with the lines that carry them.

| item | status | where |
| --- | --- | --- |
| radiation force on the gas, per x1 face, half to each cell | yes | `rad_m1_implicit.cpp:8463-8513` (dm1/dm2/dm3), write-back `:8589-8596` |
| energy exchange from the assembled row (e_gas + E conserved) | yes | `rad_m1_implicit.cpp:8448-8458` (`:8455`) |
| gas work W = vbar.dm, removed from E | yes | `rad_m1_implicit.cpp:8514-8531` |
| F in physical (r, theta, phi) components = hydro momentum basis | yes | `rad_m1_sph.cpp:26-31` |
| comoving cell flux F = F0 + (v + v.D) E with the stage velocity (vimp) | yes | `rad_m1_implicit.cpp:8538-8573`; sp rows `ImplicitVimpBuild :5756` |
| enthalpy flux (implicit_enthalpy = plm) with A_f/V_i on sp | yes; not exact on a stretched r grid | `rad_m1_curvilinear_design.md` sect. 8 |
| etotgrav (rho Phi inside u(IEN)) | yes | `rad_m1_opacity.cpp:89`, `rad_m1_implicit.cpp:8446` |
| force_reference = wb_arad (residual momentum, full work) | yes | `rad_m1_coupling.cpp:38-58`, `rad_m1_implicit.cpp:8482` (dmref) |
| a generic point-mass source | **no** | `srcterms.cpp:40-44` has const_accel only; pgens supply gravity |
| sp_cart_polar_momentum / sp_cart_all_momentum with M1 | **no** (never read by rad_m1) | `grep sp_cart src/rad_m1` empty |
| a user callback for the implicit x1 face BC | **no** | `SetImplicitX1BC` (`rad_m1_implicit.cpp:5501`) sets a type + flux only |
| the momentum deposit on sp | plain 1/2-1/2 face split, no area weight: O(dr^2) vs the cell force | `rad_m1_implicit.cpp:8463-8480` |

So the coupling needed nothing; the gap was gravity, an IC, walls and the table hand-over, which
are problem-generator work.

## 2. Route: (b), as a new file

`src/pgen/tests/rad_m1_wedge.cpp`, `<problem>/m1_test = sph_wedge` (pgen `rad_m1_beam`, the
PROBLEM-less binary), dispatched from `RadiationM1Tests2`, declared in `pgen.hpp`. Reasons:
red_giant.cpp (4900 lines) is built around two-stream/conduction/MLT with active defaults (sponge on,
required opac_table/teff/ptop, its own IC march), so an M1 mode would have to neutralise dozens of
switches; the M1 test pgen already owns every M1 hook (force reference, face-flux IC, M1 ghosts, T-S4),
and what it lacked is small and copied from proven code (RedGiantGravity's WB branch,
box_convection's table reader).

What the new setup does (all keys new, `wg_*`; nothing else changes):
- point mass `wg_gm`; gravity source = RedGiantGravity's form (plain `-rho g` + work, or with
  `wellbalance_dynamic + wb_x1` the background's area-weighted pressure drop);
- `wg_ic = grey`: grey (opacity const) atmosphere in hydrostatic + radiative equilibrium, RK4 from
  the top face inward, Eddington closure, `E_top = F_top/(c q)`, ideal gas (test A);
  `wg_ic = file`: `r rho eint [E]` column (any gas-only EOS);
- `wg_phi_eff` (default = force_reference == wb_arad): the x1 WB pair gets
  `Phi_eff = Phi - int a_ref dr`, `a_ref = kappa_t F/c` of the IC, and the M1 coupling applies only
  the residual (`SetForceReference`, face form) -- the spherical analogue of box_convection's
  `wb_phi_eff + force_reference = wb_arad`;
- `<rad_m1>/opacity = table`: `wg_opac_table` + `wg_planck_table` handed over with
  `SetOpacityTables` and the units cross-check (box_convection's code);
- closed x1 walls (`ix1_bc = ox1_bc = user`): ghosts = IC profile scaled by the adjacent active
  cell's ratio (rho, eint), v_r mirrored; M1 ghosts copy (inner) / dark (outer); the implicit
  solve uses its own face BCs (flux in, Marshak out);
- `wg_seed` (deterministic eint seed on global indices), `wg_spot_*` (isochoric hot spot);
- user history: comoving face luminosities (first face, middle, tau = wg_tau_int, Marshak face),
  L_in, E_rad, e_gas and interior (tau >= wg_tau_int) mass, KE, radial KE and momentum, sum p dV.

CUDA-safe: no lambdas inside kernels, host-built columns deep-copied to the device, no host reads of
device Views, no class members inside kernels (all captured by value first).

RESULTS: sections 3-6 (filled in below).
