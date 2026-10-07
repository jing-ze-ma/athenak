# Plaskett's star progenitor: stream-gainer setups (2026-10-07)

Source of the parameters: CANDIDATES.md section 7.2, Wade+2026 (arXiv:2609.25526).
Pgen: ry_per_accretor envelope mode, branch accretor-1006, commit ea04cd19.
Geometry and grid numbers come from the inline computation recorded in this note.
Every density ratio and temperature marked ESTIMATE is mine.

## Epoch: the start of accretion (Case A)

I chose this epoch because it is the only state on the track that the paper quotes. The quote is: "when the
secondary starts accreting it has M≃16M⊙ and R≃9R⊙. The accretion rate is ˙M≈10−4M⊙/yr" (l.710-712).

Initial binary: "M1,i = 18.2M⊙ M2,i = 16.8M⊙ Pi = 3.69 [d]". I take P = 3.69 d.

| parameter | value | status |
|---|---|---|
| gainer mass | 16 Msun | quoted |
| donor mass | 18 Msun | ESTIMATE (18.2 less winds) |
| a | 32.56 Rsun | Kepler |
| d_L1 | 15.88 Rsun | computed |
| R_acc | 9.0 Rsun | quoted |
| R_acc/a | 0.28 | |
| r_min/R_acc | ~0.17 | table, direct impact |
| GM/R | 3.39e5 (km/s)^2 (v_K 582 km/s) | |
| g | 3.77e4 code (5.4e3 cgs) | |
| P_orb | 0.458 code | |

Stream:
- Donor Teff 33 kK (ESTIMATE), mu 0.62, so c_s,don = 21 km/s. The launch speed and width come from the pgen's
  ballistic integrator.

Density ratio (ESTIMATE, Mdot 1e-4):
- Stream density at impact: Mdot/(v A), with A ~ (c_s,don/Omega)^2 = (1.53 Rsun)^2 and v ~ 650 km/s, gives
  ~9e-9 g/cc.
- Photosphere: P_ph ~ (2/3) g/kappa_e = 1.1e4 cgs at c_ph,real ~ 20 km/s, i.e. rho_ph ~ 2.6e-9 g/cc.
- So env_rho_ph = 0.3.

## Option 1: hot-surface shortcut (input bin/plaskett_env9.athinput; NOT YET RUN, submit blocked)

Grid:
- Log grid r = 6.304 (0.70 R) .. 13.5015 (0.85 d_L1), dln r = 3.10e-3.
- 246 x 4 x 2048, i.e. 2.0e6 cells.
- R_acc = face 115; dr(R_acc) = 0.0279 = r dphi.

Temperatures:
- Hot surface env_cs_ph = 64.8 km/s, so H_p = 4 dr. The real value is ~20 km/s, with H_p = 0.0106 Rsun.
- Stream env_cs_stream = 21 km/s.
- env_r_spin = R_acc; spin 1.
- Ambient as before: c_amb 300 km/s, rho 1e-6 at R_acc.

Expected cost (ESTIMATE):
- dt ~ 0.3 x 0.0279/800 = 1.0e-5 code, so ~4.4e4 cycles per orbit.
- At 2.2e8 zone-cycles/s per node, ~7 min per orbit on 1 apu node. Half an orbit takes ~4 min (apudev).

Physical compromise (ESTIMATE; stream ram pressure at ~800 km/s against the photospheric pressure):

| quantity | real surface (c_ph 20 km/s) | hot surface (c_ph 64.8 km/s) |
|---|---|---|
| ram/surface pressure | 1*800^2/(0.3*400) = 5.3e3 | 1*800^2/(0.3*4200) = 510 |
| penetration depth (n = 3 envelope, P = rho_ph c_ph^2 (c^2/c_ph^2)^4) | ~0.32 Rsun (psi 1.2e4) | ~1.7 Rsun (psi 6.3e4), i.e. 0.19 R_acc |
| envelope density contrast to r_in | 1.6e5 | 900 |

The ram pressure dominates in both cases. But the hot envelope's pressure rises much more slowly with depth,
so the stream is stopped ~5x deeper (0.19 R_acc instead of 0.04 R_acc). This is the main physical compromise
of option 1 for Plaskett: the AM is deposited deeper and into a lighter envelope.

## Option 2: resolved real surface (hydro), design and cost (ESTIMATE)

Real surface: c_ph 20 km/s, H_p = 0.0106 Rsun; dr_min = H_p/4 = 2.65e-3 Rsun.

The fine zone is only needed over the top ~0.1 Rsun of the star and the atmosphere, R_acc - 0.1 .. + 0.1. Below
it the n = 3 envelope's own H_p ~ H_ph + depth/4 is already 3-4 log cells per H_p at depth 0.3 Rsun.

Grid (grid_design.py approach, polynomial + plateau, no core change):
- ~75 fine cells plus 2 x ~25 transition cells (ratio <= 1.1) on top of the 246-cell log grid.
- nx1 ~ 360; with nphi 2048, 360 x 4 x 2048 = 2.95e6 cells.

Time step and cost:
- dt = 0.3 x 2.65e-3/(800 + 20) = 9.7e-7 code, i.e. 4.7e5 cycles per orbit.
- At ~2.3e8 zone-cycles/s per node: **~1.7 h per orbit on 1 apu node, ~0.9 h on 2, ~0.5-0.6 h on 4.** At this
  size, 8 GPUs scale below ideal.
- This is ~9x cheaper than the resolved RY Per photosphere (15.5 h/orbit/node), because H_p/R is 1.2e-3 here
  against 8e-4, the cells are fewer, and P_orb is 3.69 d against 6.86 d.

## Radiation (M1, he_star_m1 style) on top of option 2

Requirements:
- **Physical units.** The run is no longer scale-free: Mdot = 1e-4 Msun/yr sets the density in g/cc. The donor
  and gainer temperatures enter through T.
- **Coupling.** <rad_m1> with the implicit transport, coupled to hydro. The etotgrav, wellbalance_dynamic and
  wb_x1 machinery that he_star_m1 requires is already in place.
- **Opacity.** Electron scattering (0.34 cm^2/g) plus a Rosseland/Planck table (OPAL/Kramers) for T > 2e4 K.
- **Radiation boundary conditions.** The stellar flux L(r_in) enters at the inner wall: today's log L 5.1 is a
  stand-in, and the gainer at 16 Msun has log L ~ 4.5 (ESTIMATE). The outer boundary is free streaming. The
  stream injects T_don.
- **Pgen work.** A he_star_m1-style radiation-energy initial state for the envelope, i.e. radiative HSE with the
  simulation EOS and opacity (user rule). Gamma_e is ~0.1 now, rising to 0.45-0.7 during rapid accretion
  (Wang+2026), so radiation pressure must enter the HSE and the force.
- **CFL.** rad_signal_speed (implicit transport, so no light-speed limit) and cfl 0.3 for radiation-dominated
  runs (user rule 10-01).

Cost (ESTIMATE; the factor is not measured for this geometry):
- The implicit M1 typically costs 3-5x the hydro per cycle on the He boxes.
- **~5-9 h per orbit on 1 node, ~1.5-2.5 h on 4 nodes.**
- Plus, on the order of a week of setup work: units, opacity, inner-flux BC, radiative HSE initial state, and
  test gates (static radiative envelope, then the stream).

## Option 2 runs (2026-10-07)

**Setup.**
- Input: bin/plaskett_env10.athinput.
  - Grid: 500 x 4 x 2048 cells.
  - Fine cells: dr 2.5-2.8e-3 Rsun on r 8.80..9.13.
  - Surface: c_ph 20 km/s; stream: 21 km/s.
  - WB everywhere, nghost 2, no fofc.
  - r_acc: 9.00129, a cell face.
- Binary: env7 (md5 a5ad0a18c1d3a73e4cb41e833f1821d5).
- r_top: 9.1815.

**Smoke (job 12109320).** Clean in both arms:

| arm | dt | zone-cycles/s per node |
|---|---|---|
| with stream | 1.94e-6 | — |
| without stream | 1.83e-6 | 2.43e8 |

**Envelope only, 0.2 orbit (job 12109350, q10_a1, 1 apu node, 13 min): PASS.**
- 49 900 cycles, dt 1.84-1.85e-6, no collapse.
- Envelope (r < R_acc):
  - max \|dv\| 0.024 km/s;
  - Menv drift -8.3e-6;
  - Jenv drift -1.6e-5.
- Atmosphere (rho > 1e-3 rho(R_acc)):
  - max \|dv\| 7.8 km/s at r 9.053;
  - max \|drho/rho\| 3 %.
- Note: scripts/ana_env.py's shell radii are hardcoded for RY Per (R_acc 4.06), so its Omega lines do not apply here.

**Stream run (job 12109351, q10_s1, 2 apu nodes): FAILED, dt collapse at t = 0.03531 (cycle 19177, 0.077 orbit).**
- Up to the last dump (t = 0.0345), dt is 1.84e-6, limited by hot ambient cells at r 9.08 (c_ad 392 km/s).
- At t = 0.0345 the stream has reached r ~ 9.2 near phi 2-13 deg.
  - In that sector, dense gas sits above R_acc: rho 47 at r 9.20, 159 at 9.10 and 255 at 9.03, with T 2.5-3.7e3 (c 50-60 km/s).
  - This is the tidally bulged envelope top on the donor side. R/R_L ~ 0.75, and the initial state follows the full tidal Phi_wb.
- The collapse falls between dumps.
- The dense restart diagnostic (qd1/qd2) was blocked by the permission check; its command is in PENDING_CMDS.txt. The mechanism is NOT yet confirmed.
- Following the user's instruction, there was no continuation, no spin-7.2 twin and no physics change.

### Collapse diagnosis (job 12109481: qd1 restart t 0.0229 -> 0.0344, then qd2 with dumps every 2e-5)

- Deterministic: qd2 collapses at the same cycle (19177, t = 0.035314) as q10_s1.
- Up to the last dump (t = 0.035312), dt is 1.82e-6. The limiting cell is hot ambient at r 9.11-9.12, phi 68 deg (rho 7e-7, v_r -60 km/s, c_ad 380 km/s). There is no precursor and no dt decline.
- The collapse happens in ONE cycle. dt drops to 1.4e-9 at r = 9.1635 (R_acc + 0.16, i.e. inside the WB zone, r_top = 9.18), phi = 61.2 deg. That cell goes to the floor density 1e-7 with \|v\| ~ 2.6e6 km/s: a positivity failure.
- The site is not the donor-side bulge (phi 2-13 deg). It is the thin layer above the surface at phi ~ 61 deg. There, stream-fed gas sweeps along the surface: rho ~1e-2, T ~ 400, v_phi ~ 330 km/s, v_r -10 to -250 km/s. The hot ambient above it has rho 5e-7, T 8e4 and is nearly at rest.
- Across this interface the density contrast is ~2e4 and the shear ~300 km/s, over 3-5 radial cells (dr 3e-3).
- The interface also shows strong theta structure inside the 4-cell band (the 4 j cells differ by 10-100x in rho and by 60 against 250 km/s in v_r), and isolated cells far below both relaxation targets (T = 4.7 next to 3e4).
- **Mechanism:** a one-step positivity loss in a cold-dense/hot-thin shear layer inside the WB zone, under a strong shear. It is the same class as RY Per's t = 0.069 failure (WB/plain face, Mach-10 shear) and 0.092 (FOFC).

### Candidate fixes (for the user's decision; none applied)

| fix | type | expected effect | caveat |
|---|---|---|---|
| hydro/fofc = true (+ mesh/nghost = 3) | input keys | catches exactly this one-step positivity loss. On RY Per it moved the collapse from 0.069 to 0.092; it may hold here, where the real surface is resolved | — |
| hydro/wb_rmax = R_acc | input key | WB stays on the star only; on RY Per it removed the impact collapse | the thin layer above the surface would then be plain PLM |
| lower the density contrast at the interface (env_amb_rho 1e-6 -> ~1e-5, P_amb 0.4 against the cold layer's ~4) | input key | gentler contrast, smaller positivity risk | it changes the ambient criteria |
| per-cell WB switch under strong shear or inflow; a positivity-preserving (density/pressure) limiter on the reconstruction | core change | targets the mechanism directly | needs a code change |
| a true r-phi 2-D run instead of the 4-cell theta band | not supported as is | would remove the theta structure | AthenaK needs nx2 > 1 whenever nx3 > 1 |

**Recommended first try:** fofc + nghost 3 (no code), combined with wb_rmax = R_acc. Each is cheap on apudev up to t = 0.04.

### fofc + nghost 3 + wb_rmax = R_acc (plaskett_env11, job 12115722)

The smoke is clean. q11_s45 reaches t = 0.045 (cycle 24177) with 0 COLL and 0 FATAL, i.e. past the old
collapse at 0.0353.

**dt.** It falls from 1.8e-6 to 9.3e-7 by t = 0.045 as the stream arrives. The dump estimate of the
limiting cell is hot ambient at r 9.09 (c_ad 395 km/s).

**FOFC firing** (t 0.0445, cumulative over the window):
- at-floor cell-stages: 8.3e7;
- other cell-stages: 8.8e5.

Almost all the firing is in cells at x1 index >= 327 (r >= 9.34), i.e. floor/ambient gas above the surface.

**Next step.** The 0.5-orbit run (q11_s1, 2 nodes) and the envelope-only check (q11_a1) are blocked by the
permission check. Their commands are step 2 in PENDING_CMDS.txt.
