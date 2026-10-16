# TASK 2026-10-16 caltech -> raven: accretor bigdom HANDOVER (Raven now owns the code and the tests)

**User decision (10-16): the whole accretor bigdom project moves to Raven.** Caltech has stopped its accretor bigdom
testing. From now on Raven OWNS the branch accretor-bigdom-1015: the code, the tests and the reports.
**Report back by NOTE. Merge decisions stay with the user.** Caltech keeps only the orbit1 production run
(branch accretor-rhd-1015, binary fdb18517) and does not touch this branch any more.

Replies to: NOTE-2026-10-16-raven-accretor-bigdom.md (incl. ADDENDUM 1), NOTE-2026-10-16-raven-accretor-3d.md,
NOTE-2026-10-16-raven-accretor-donorlight.md (3bd0da2f), NOTE-2026-10-16-caltech-accretor-fix3-yes.md (32540c1a).
This replaces the planned "v4" TASK.

## Branch, sha, build
- Branch **accretor-bigdom-1015 @ 99adf52d** (code **5b0d4bc1**; 99adf52d adds docs only). Design note:
  docs/dev/accretor_rhd_design.md **section 15** (v4), sections 13-14 (bigdom, 3-D).
- v4 code commits: 5ebfac82 (all four features) + 5b0d4bc1 (sink_zero_v kept out of the parameter dump).
- Build: your recipe (gcc/13 + cuda/12.6 + openmpi_gpu/5.0, Kokkos 4.6.02 CUDA AMPERE80, Release, MPI,
  `-D PROBLEM=ry_per_accretor -D CMAKE_CXX_FLAGS=-ffp-contract=off`). **Full fresh build**: src/pgen/pgen.hpp
  changed (new hook pointer), and src/rad_m1/rad_m1_coupling.cpp calls it.
- It builds: the Caltech GPU nofma build of 5b0d4bc1 is BUILD_OK (H200, md5 064f2963). **Nothing of v4 has run.**
  The Caltech static check and gate were queued but stopped before they started (handover).

## Deck and tables
docs/handover/accretor-bigdom-1015/ on the branch (unchanged): `bigdom.athinput` (v3, md5 3a818e64),
`column_z0.02_p381.txt` (78f4a632), `rosseland_ext2_gs98_x0.70_z0.02.txt` (ad516c7e),
`planck_ext2_gs98_x0.70_z0.02.txt` (99d668f6), and `b3d_half.athinput` (3-D). `sed "s#@HDIR@#$PWD#"` as before.
Command-line overrides only work for keys already in the deck, so put new keys (`dt_report = true`,
`donor_irr_scale = ...`, `donor_gdark = ...`) into the <problem> block of the deck file itself.

## State of the user decisions
| item | state |
|---|---|
| Pinned-rate nozzle at r 14.5 (stream_mdot) | DONE (v3, tested by you) |
| Donor L via donor_lum (69183 Lsun) | DONE |
| 2-D donor light scaled to the 3-D dilution (problem/donor_irr_scale) | IMPLEMENTED, UNTESTED, NOT CALIBRATED (below) |
| Von Zeipel gravity darkening of the lobe (problem/donor_gdark, default 0.25) | IMPLEMENTED, UNTESTED |
| Back-sink sector holds the STREAM's own light (nozzle T profile), not donor light | IMPLEMENTED, UNTESTED |
| Sink-cell momentum zeroed after the radiation step (problem/sink_zero_v, default true) | IMPLEMENTED, UNTESTED |
| Fix 3: vet_col lateral-beam test | APPROVED by the user (32540c1a); yours to run |
| Fix 4: force/sponge guard for atmosphere gas | NOT decided |
| 3-D physics (vertical gravity, Coriolis/centrifugal in 3-D, IC, WB reference, ambient; theta walls) | TO DO; you said you would wait for it before any longer 3-D run |

### Donor light: which form
The implementation is the **scaled lobe** (your stopgap form), on top of the gravity-darkened lobe, as our user
decided (scale + gravity darkening + stream-lit sector). Your preferred physical form (an inward flux
L_d/(4 pi d^2) only within the donor's solid angle) is NOT implemented; as the owner you may choose it instead.
- `donor_irr_scale` multiplies the held lobe E (uniform or gravity-darkened). Given: used as is. Absent: **1 in
  3-D**, **auto in the 2-D band**: s = s_geo x s_m1.
  - s_geo = int T_loc^4 dA / (2 D oint T_loc^4 dl): the band with reflecting theta walls is a donor infinitely long
    in z, whose light dilutes as 1/D (cylinder) instead of 1/D^2 (sphere). For a sphere s_geo = R/D.
    D = a - donor_irr_rcmp (default 9.6) = 23.66. **s_geo = 0.636 with beta 0.25** (0.523 uniform). The start-up log
    prints s, s_geo and the 3-D target flux L_d/(4 pi D^2) = **1.52e5 code** (1.54e5 uniform).
  - s_m1 = `donor_irr_m1fac` (default 1) is the band M1/vet_col excess over ideal cylindrical dilution. **It is
    not calibrated.** Your ~30x at s = 1 is ~15x more than ideal transport would give (1/s_geo ~2x), so expect
    s_m1 ~1/15 and s ~0.035-0.04, close to your 1/30 stopgap (T_d,eff ~11 kK). The default auto (s_m1 = 1) is
    therefore **too bright by ~15x** until calibrated: set donor_irr_scale explicitly (e.g. 0.033) or calibrate.
  - Calibration: two static arms to t 0.00946 (donor_irr = false; and donor_irr_scale = 1): facing-side donor flux
    at r 9.6 (|phi| < 10 deg) = -(F_r(irr) - F_r(no irr)); then s = 1.52e5 / that flux. M1 transport is linear in
    the held E, up to the gas coupling. Then set donor_irr_m1fac = s / s_geo, so that auto gives s.
  - Part of the excess may be the closure (vet_col builds the Eddington tensor from radial columns, so the lateral,
    azimuthal donor light is not represented): fix 3, your lateral-beam test.

### Gravity darkening (donor_gdark = beta, default 0.25; donor_gdark_tfloor, default 0.3; 0 = old uniform T_d)
- Surface: the L1 equipotential of the rotating-frame Roche potential the code uses. g_eff = |grad Phi|.
  T_loc^4 = C g_eff^(4 beta), floored at (0.3 T_mean)^4. Each lobe cell holds E = s a T_loc^4, F = 0, with g_eff at
  the surface point where the ray from the donor centre through the cell leaves the lobe (rays through the neck:
  g = 0 -> the floor).
- Normalisation: int sigma T_loc^4 dA = L_d over the **full 3-D lobe surface**, also in the 2-D band (the band
  holds the equatorial cut of the 3-D donor). A_lobe 2133 Rsun^2, T_mean 25.9 kK.
- Expected in-plane T_loc by angle from the L1 axis (python twin): 7.8 kK at the neck (floor), 16.9 kK at 5 deg,
  21.4 kK at 15, 24.2 kK at 30, 26.2 kK at 60-90 (side), 25.5 kK at 120, 24.2 kK at 180 (back, outside the
  domain). Check the start-up table against these.

### Back-sink sector (r above the nozzle, |dphi| < 3 sigma + 1 sigma)
Old: held a T_d^4 (doubled Eirr in the stream phase). New: within the nozzle width (|dphi| < 3 sigma; in 3-D also
|z_f| < 3 sigma_z) it holds E = a T^4 with the nozzle's T profile (core env_t_stream, diffusion profile, floor
stream_t_surf ~19 kK), F = 0. Outside it its radiation is left to the transport (dark, like the ambient). Only with
stream_rad = true. Not scaled (stream light, not donor light). Booked with the nozzle's radiation (Ebnd). The
one-zone T(r) back toward L1 is not in the code, so the nozzle profile is held unchanged.

### Sink momentum (sink_zero_v)
New hook user_m1post_func, called by RadiationM1::HydroConToPrim after the M1 coupling and before the fluid
refresh. It zeroes the momentum of every sink cell (lobe and back sector, not the nozzle), removes the KE from E_tot
and books it in Emsk. Active only with donor_irr on or a live back sink (stream on). The hydro dt is computed before
the M1 stage, so dt should be unchanged. It should remove your r 14.6 row artefact.

## Tests to run (yours; test only, report by NOTE)
1. **Regression gate first (rule: before any merge, the old 2-D decks must stay bitwise unchanged):** 5b0d4bc1 vs
   08625ed7 (or your efe33e6a binary), 2-D v3 deck + dt_report, donor_irr = false, new keys absent, 20 cycles:
   hst, bins (hydro_w, m1, including the PAR_DUMP header) and cycle/dt lines identical. Also run your earlier
   old-deck gate set if you can.
2. Startup check: the printed gdark table, s_geo and the target flux match the numbers above.
3. Calibration of the scale (above), unless you take your directed-flux form.
4. **2-D v4 through impact (tlim 0.036)**: the v3 deck with the calibrated scale, the static phase to 0.00946, then
   the stream on by restart. Pass checks:
   - no FATAL; no NON-CONVERGED storms; list the dt retries;
   - pinned rate: Mnzf rate (col 25) = 19.795 within 0.5 % after the switch-on; Min - Mnzf constant after the fill;
   - **dt >= ~6e-6 before impact** (reference 7.03e-6; v3 fell to 2.4-3.3e-6);
   - thin-atmosphere (r 9.1-9.6) max |v| at t 0.005-0.0095 against your no-irr 47 and irr 228 km/s, and through
     the stream phase against v3 (630-940 km/s);
   - facing-side donor flux at r 9.6 against the 1.52e5 target (within ~1.5x);
   - mask and sinks clean: lobe at floor, the r 14.6 back-sink row at rest in the outputs, Mmsk/Emsk smooth;
     Eirr rate in the static phase (back sink off) and in the stream phase (the back sector no longer adds donor
     light);
   - arrival phi and time of the first impact against orbit1 (15-16.5 deg).
5. Fix 3 (vet_col lateral beam), as approved.

Your v3 impact job (31212975) runs the old over-irradiated deck. Stop it or let it finish for comparison, your choice.

## 3-D
Keep the scale at 1 in 3-D (the default) unless the 3-D benchmark shows the same over-irradiation. The 3-D physics
items above are still to do.
