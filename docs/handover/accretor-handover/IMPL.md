# RY Per accretor: implementation and tests (2026-10-06)

## What was built

The code is on branch `accretor-1006`, off rt-integration 495df8fc, in the worktree
/viper/ptmp2/jinma/wt_accretor. It is committed but not pushed or merged.

| commit | content |
|---|---|
| b226528c | **core fix**: isothermal pressure in the sp / cubed-sphere geometric sources (coordinates.cpp, 5 sites) |
| f5dfc1ef | new pgen `src/pgen/ry_per_accretor.cpp` + `inputs/hydro/ry_per_accretor.athinput` |
| 1e8d0ff8 | **core fix**: Driver periodic NaN scan read u(IEN) for isothermal (out of bounds) |

Two pre-existing bugs blocked an isothermal spherical-polar run. Both fixes leave the ideal and general EOS paths bitwise unchanged:

1. `SrcTermsSphericalPolarHydro` (and the MHD, WB, gnomonic and cart-rows variants) computed p = (gamma-1) w0(IEN). An isothermal EOS has no IEN slot, so the 2p/r term was garbage. A hydrostatic isothermal atmosphere then fell in at 2 c_s^2/r (measured: v_r grew uniformly, a ~ 2 c_s^2/r). The fix uses p = rho c_s^2 when `!is_ideal`.
2. The NaN scan in `Driver::Execute` tested u(m,IEN,...) unconditionally.
   - On GPU this was hipErrorIllegalAddress at cycle ~90-100 of every isothermal run (smoke job 12105143; diag job 12105223, all 5 variants).
   - The CPU Debug build gives: `cons [0,4,0,0,0] extents [1,4,...]`.
   - The fix tests the density when there is no energy slot.

Builds: own build dirs under /viper/ptmp2/jinma/accretor_1006, made by `build.sh cpu|gpu72` and `build_dbg.sh`.
- GPU: gcc/16 rocm/7.2 hipcc 7.2.4, GFX942_APU, MALLOC_ASYNC off.
- Binaries in `bin/`:
  - `athena_ryper_gpu72_1e8d0ff8` (md5 b6ab7220c3dab29fb93bfd03d5bfc6b5);
  - `athena_ryper_cpu_1e8d0ff8`;
  - the matching input `ry_per_accretor_1e8d0ff8.athinput`.
- Style: `cpplint` is clean on ry_per_accretor.cpp. driver.cpp and coordinates.cpp have no new findings (the pre-existing count is unchanged). There are no tabs, no trailing whitespace, and the files are mode 644.

## Physics, units, keys

Units:
- length Rsun;
- velocity km/s (time unit 6.957e5 s, P_orb = 0.8519 code);
- density in units of the stream peak density (scale-free).

Frame and potential:
- Omega = sqrt(G(Ma+Md)/a^3), from Kepler. That gives 1.0601e-5 /s, against 2 pi/P = 1.0595e-5 (0.05 %). With Kepler, L1 is exact.
- Phi = -GMa/r - GMd/|r - a e_x| - Omega^2 |r - x_cm e_x|^2 / 2, with x_cm = a Md/(Ma+Md) = 6.457 Rsun.
- **Indirect term:** expanding the CM centrifugal term about the accretor gives a constant -Omega^2 x_cm e_x = -(G Md/a^2) e_x. That is exactly minus the accretor's acceleration toward the donor, so no separate indirect term is needed. This is derived in the file header.
- Gravity is the face-difference gradient of Phi at theta = pi/2 (in-plane, no theta component), precomputed per (m,k,i).
- Coriolis is explicit, in-plane: a_r += 2 Omega v_phi, a_phi -= 2 Omega v_r.

EOS and the reused energy machinery:
- Isothermal, `problem/thermo = isothermal`. `adiabatic` is reserved and FATAL; what it needs is in the header: etotgrav phi0 fill like he_star_m1, and the 09-28 rot_potential double-count guard.
- With no energy equation, the etotgrav and rot_potential machinery does not apply. Only the Coriolis pattern of dhj_rt and the HeStarBC skeleton are reused.

Boundaries:
- **Inner** (x1min = R_acc): absorbing diode, v_r = min(v_r, 0).
  - Density: `inner_rho = copy` (default) or `hse`.
  - Spin: `inner_slip = noslip` gives v_phi = (spin-1) Omega r; `free` copies v_phi.
  - `inner_vr = wall` is a closed mirror wall, for tests only.
- **Outer** (x1max = 0.85 d_L1 = 16.241 Rsun):
  - Stream window |phi - phi_s| < 3 sigma, with a Gaussian rho and the ballistic v_r, v_phi. These come from a C++ RK4 copy of roche_stream.py's integrator (same launch from L1 at c_s,don), run at start-up.
  - Elsewhere: v_r = max(v_r, 0), and density from `outer_rho = hse` (default) or `copy`.
- Theta: reflecting band of 4 cells (dtheta = dphi). Phi: periodic.

Stream state used at r_out (printed at start-up):

| quantity | value |
|---|---|
| phi_s | 4.157 deg |
| v_r | -72.58 km/s |
| v_phi | +35.44 km/s |
| width c_s/Omega | 0.8643 Rsun |
| sigma_phi | 0.0532 rad |

The C++ integrator reproduces roche_stream.py:
- r_min = 2.6946 Rsun (script: 2.69);
- at R_acc: phi 73.74 deg, v_r -350.4, v_phi 545.1 km/s, 57.3 deg from the normal (script, sampled at r = 4.02: 74.6 / -349 / 551 / 57.6).

Problem keys, all with defaults in the input: m_acc 6.24, m_don 1.69, a_sep 30.3, period, r_acc 4.06, t_don 6250, mu_don 1.27, spin, inner_slip, inner_rho, inner_vr, outer_rho, hse_cap, stream, rho_stream 1, rho_amb 1e-5, stream_nsig 3, stream_width, stream_phi_deg/stream_vr/stream_vphi (overrides), init ambient|hse, thermo.

History (`ryper.user.hst`, code units, summed over ranks):

| column | meaning |
|---|---|
| Mdom | mass in the domain |
| Jdom | inertial Lz in the domain |
| Min | cumulative mass in through the stream window |
| Mout | cumulative mass out through the rest of r_out |
| Macc | cumulative mass accreted through r_in |
| Jacc | cumulative inertial Lz into the star |
| Jstr | the part of Jacc not carried by the edge-cell v_phi (the no-slip wall stress; approximate, uses the cell-centred v_phi) |
| Jout | cumulative inertial Lz out through r_out |
| d* | the same six as rates over the last history interval |

- Cumulatives are time integrals of the Riemann fluxes with the RK stage weights; only rk1/2/3 are supported.
- They restart from 0 on a restart.

## Tests

### (a1) At-rest atmosphere

No stream; synchronous, no-slip star; `init = hse`, hydrostatic inner and outer ghosts; 0.1 orbit (t = 0.08519). CPU, quarter resolution 112x4x512, 8 ranks.

- **Closed wall** (`inner_vr = wall`):

  | c_s (km/s) | max\|v\|, whole domain | max\|v\|, inner 5 cells | where the maximum sits |
  |---|---|---|---|
  | 300 | 2.7e-4 c_s | 1.7e-5 c_s | — |
  | 150 | 3.0e-2 c_s | 1.7e-4 c_s | at the open outer boundary, r 15.5, low density |

  rho_max is unchanged (1.1141e-1 vs 1.1131e-1 for c_s 150). Directory: tests/a1w2_hot{300,150}.
- **Absorbing diode** (production inner BC), c_s 300: max\|v\| = 3.5e-2 c_s. At c_s 150 the atmosphere drains through the absorbing surface: Macc 0.29 of 0.30 by 0.1 orbit, accelerating. That is expected, because an absorbing surface gives no support. Directory: tests/a1_hot{300,150}.
- The test with the production c_s = 15.5 km/s is not possible: H/R ~ 1e-3 is unresolved, so the ambient must fall in.

### (a2) Cold ballistic stream

c_s = 1 km/s, GPU. Jobs 12105359 (r_in 4.06, full resolution 448x4x2048) and 12105358 (r_in 2.0, 336x4x1024), both 0.3 orbit, analysed with tests/ana.py.

With r_in = 4.06:
- The ridge phi(r) follows the ballistic orbit within -0.42 to +0.45 deg from r = 14.9 down to the surface.
- Impact (peak inward mass flux), run vs ballistic:

  | quantity | run | ballistic |
  |---|---|---|
  | phi | 74.09 deg | 73.65 deg |
  | v_r | -349.3 km/s | -350.4 km/s |
  | v_phi | 545.7 km/s | 545.1 km/s |
  | angle from the normal | 57.4 deg | 57.3 deg |

- 100 % of the inward flux lies within ±10 deg.

With r_in = 2.0:
- The ridge is within 0.5 deg of ballistic down to r = 4.4, and 0.8 deg at r = 3.86.
- The minimum r of dense gas (rho > 0.3) is **2.774 Rsun vs ballistic r_min 2.695 (+3 %)** at t = 0.10-0.15. The cell width there is 0.017.
- Near periastron the radius-based ridge tracker is ill defined, and the stream later self-intersects.

### (c) Spin 1.0 vs 7.2

Production c_s 15.5, stream on, CPU quarter resolution, 0.5 orbit (tests/c_spin1, c_spin72). Means are over t = 0.25-0.426 (tests/budget.py):

| arm | Mdot_in | Mdot_acc | j_acc / j_Kep(R_acc) | Jstr/Jacc | dJstr (last interval) |
|---|---|---|---|---|---|
| spin 1.0 | 125.0 | 125.1 | 1.0647 | 4.0e-4 | 110 |
| spin 7.2 | 125.0 | 125.1 | 1.0647 | 3.0e-4 | 90.6 |

- The sign is sensible: the faster star takes about 18 % less wall-stress torque.
- The size is negligible. The accreted angular momentum is the ballistic stream's own (1.065 j_Kep vs ballistic 1.067). The wall stress is about 3e-4 of the advected torque, because the impact is supersonic and the diode surface absorbs it.
- **So with an absorbing rigid surface the star's spin hardly affects the torque.** A spin dependence needs either a resolved envelope or boundary layer (DESIGN stage 2) or a viscous or non-absorbing surface.
- These runs used binary 1e58352a, which predates the driver fix. The fix only touches the NaN scan.

### (b) GPU smoke

apudev, 2 x MI300A, production input, nlim 300, job 12105298: rc 0, 0 FATAL, 0 NaN.

Throughput, measured:

| run | zone-cycles/s per node | per GPU |
|---|---|---|
| smoke | 3.41e8 | 1.71e8 |
| a2_rin406, 37677 cycles | 3.44e8 | 1.72e8 |

The assumed 3e8-1e9 per GPU is not reached. The likely cause is the 4-cell theta band, which carries 2 ghost cells on each side.

Cost:
- Steady dt at full resolution with the stream on the star: 6.58e-6 code (4.6 s), so 1.29e5 cycles per orbit.
- **23 min per orbit per apu node**, measured from the a2_rin406 wall time: 382.6 s for 35800 cycles.
- 10 orbits are about 3.8 h on 1 node.

A CPU Debug build (Kokkos bounds checks) at production resolution showed no out-of-bounds access in 210 cycles, through 2 NaN scans plus bin and hst outputs (tests/dbg_full2).

## Not done / open

- gamma = 5/3 is not built; the key is reserved.
- No production run was launched.
- The theta band is not well balanced against the sp cot(theta) terms. The tests show this is negligible (max\|v_theta\| < 2e-3 km/s).
- On a restart the cumulative flux integrals start again from 0 (they are not stored in the restart).
- The spin result above implies the stage-1 absorbing BC cannot measure spin-dependent accretion torque. This needs a user decision on the stage-2 inner boundary.
