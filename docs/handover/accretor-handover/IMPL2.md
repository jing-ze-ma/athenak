# RY Per accretor, stage 2: a resolved stellar envelope instead of the absorbing surface

Branch `accretor-1006` (worktree /viper/ptmp2/jinma/wt_accretor), new commits on top of 1e8d0ff8.
Numbers in section 1 come from `scripts/env_tradeoff.py`, scaled from the stage-1 measurements in
IMPL.md (dt 6.58e-6 code at full resolution, 1.29e5 cycles/orbit, 3.44e8 zone-cycles/s per apu node,
23 min/orbit). Code units: Rsun, km/s, time 6.957e5 s, P_orb = 0.8519.

## 1. Trade-off (written before coding)

Accretor: GM/R = 2.93e5 (km/s)^2 (v_K = 541 km/s), g = 7.22e4 code (1.04e4 cgs).
Run temperature: c_iso = 15.5 km/s, so H_p = c^2/g = 3.33e-3 Rsun (H/R = 8.2e-4).
Stage-1 log grid: dln r = 3.094e-3, dr(R_acc) = 0.01256 Rsun = 3.8 H_p.

| option | what it is | cells per H_p | dt vs stage 1 | cells | cost per orbit (1 apu node) | physics |
|---|---|---|---|---|---|---|
| A. genuine isothermal skin | c = 15.5 everywhere, dr = H_p/4 = 8.3e-4 near R_acc | 4 | 1/15.1 (impact \|v\|+c = 680 km/s over 8.3e-4) | ~580 x 4 x 2048 | ~15.1 x 1.3 x 23 min = **7.5 h** | skin only ~20-25 H_p = 0.08 Rsun deep: 0.7 R_acc would need e^523 in density. Its AM capacity is tiny, and the wall under it is the real star again. |
| B. **gamma 5/3 + n = 3 polytrope + relaxation (CHOSEN)** | ideal gas; envelope P/rho = c_ph^2 + psi/4 (psi = depth in the effective potential) down to r_in = 0.70 R_acc; Newtonian relaxation of T toward that profile, toward c_ph = 15.5 outside | 3.8 at the photosphere (WB-isothermal held); >= 4 below 0.19 Rsun depth (4.6 % R), ~30 at 0.7 R | 1.0 (deep: c_ad 230-240 km/s + v_phi 130 over dr 8.8e-3 gives 7.1-11.5e-6 > 6.58e-6) | 564 x 4 x 2048 (+26 %) | 1.26 x 23 min x ideal/WB overhead: **~30-40 min (estimate; measured in test b)** | real stellar-like interior (n = 3 = Eddington standard model, stably stratified for adiabatic perturbations); outside the star the gas is relaxed to the stage-1 isothermal temperature |
| C. globally hot isothermal | one c_s large enough to resolve H_p | c = 60 km/s for H_p = 4 dr; c >= 80 km/s for an e^20 envelope to 0.7 R | ~1/1.1 | 564 | ~1.4 x 23 min | changes the stream: Mach 42 -> 8, a thick pressure cushion. Rejected. |
| B' | B with t_relax -> 0 | as B | as B | as B | as B | locally isothermal T(r) envelope (neutral buoyancy); same code, one key |

Envelope depth vs dt (n = 3, spin 7.2, rotation in phi): r_in = 0.8 R gives dt 9.0e-6; 0.7 R gives 7.1e-6;
0.6 R gives 5.5e-6, which would limit the global dt by 16 %. **r_in = 0.70 R_acc**. Exactly 116 log
cells below R_acc gives r_in = 2.835492977 and nx1 = 564, so r = R_acc is a cell face (face index 116).

Density contrast across the envelope (n = 3): rho(0.7 R)/rho_ph = 2.3e6 (spin 1) and 3.0e6 (spin 7.2).
n = 1.5 gives 3e3, but it is only marginally stable (adiabatic), so n = 3 is the choice.

### What B changes physically (USER DECISION flagged)

1. The envelope is a model, not RY Per's real envelope. It is an n = 3 polytrope anchored at the stage-1
   photospheric temperature (c_iso 15.5 km/s), with surface density `env_rho_ph` in stream units.
   - ESTIMATE: P_ph ~ (2/3) g/kappa with kappa ~ 0.5 gives rho_ph ~ 6e-9 g/cc.
   - The stream density is 2.7e-10 x (Mdot/1e-6 Msun/yr) (DESIGN.md).
   - So `env_rho_ph = 20` corresponds to Mdot ~ 1e-6. This is the default.
   - The run stays scale-free only up to this single ratio, which now sets how deep the stream penetrates.
2. The run is no longer isothermal. Gas outside the envelope is relaxed to c_ph = 15.5 km/s with time
   `t_relax`; the default 1e-4 code is ~15 steps, i.e. near-isothermal like stage 1. Envelope gas is relaxed
   toward its own initial T profile with `t_relax_env`.
   - With a short t_relax_env (B') the envelope loses its buoyancy.
   - With a long or infinite one, the impact heating stays in the envelope.
   - The default is 1e-4, the same as outside, so impacts are effectively isothermal everywhere.
3. The photosphere itself (the top ~0.19 Rsun) is NOT resolved: H_p < 4 dr there. The well-balanced
   scheme holds it in equilibrium. For ram-pressure balance (rho_s v_r^2 = P_env), the stream stops only
   ~0.03 Rsun (2-3 cells) below the photosphere at env_rho_ph = 20. So the AM is deposited in an
   under-resolved skin and carried inward by numerical shear viscosity only. A resolution study of that
   layer would need option A's cost.
4. r_in wall: closed, mirror v_r; v_phi free-slip (`env_wall_slip = free`, default; wall torque ~ 0) or
   held at the initial spin (`noslip`, a rigid core).

## 2. What was built (option B)

Pgen `src/pgen/ry_per_accretor.cpp`, new mode `problem/inner = envelope`. The key is read only when
named, so a stage-1 input's parameter dump and run are unchanged: the surface kernels are untouched,
and only the hist loop bound became the explicit constant kNaccS = 6. Input
`inputs/hydro/ry_per_accretor_env.athinput`.

Requirements (FATAL otherwise):
- `thermo = adiabatic`;
- `<hydro>` eos = ideal (gamma 5/3), hllc, etotgrav, wellbalance_dynamic, wb_x1, wb_option = isothermal,
  wb_rho = true.

Grid:
- r_in = 2.835492977 (0.698 R_acc) to r_out = 16.241;
- 564 x 4 x 2048, with R_acc = log-grid face 116;
- every MeshBlock must span R_acc (meshblock/nx1 = mesh/nx1).

### Potentials

| potential | where it is used | definition |
|---|---|---|
| TRUE Phi | phicc0 / phi0 (etotgrav energy) | Roche + orbital centrifugal, as stage 1 |
| Phi_wb | the x1 WB pair (EnableWBEffectivePotential) | Phi - (s^2 - 1) Omega^2 r^2 / 2 below r_top; FLAT in r above r_top |

Phi_wb is the effective potential of a solid-body star at spin s in the rotating frame: orbital
centrifugal + Coriolis + v_phi^2/r give s^2 Omega^2 r. Above r_top it is flat: the WB background there
equals the cell state, which gives plain PLM and no WB force.

r_top is the first face above R_acc beyond `env_top_hp` = 15 photospheric scale heights (independent of
the floor; 10-07): 4.1106 at half resolution. It can be overridden with `env_r_top`.

### Forces

- **x1**: the WB pressure form (he_star_m1 pattern) + rho (-dD/dr), with D = Phi - Phi_wb as a face
  difference.
  - This is exactly -rho dPhi/dr in total, and exactly balanced for the rotating rest state.
  - Above r_top it is plain gravity.
- **phi**: tidal gravity, explicit (as stage 1).
- **Coriolis**: explicit (as stage 1).
- **No explicit centrifugal work**: the potential's work is in the etotgrav flux (09-28 guard).

### Thermal relaxation

- e -> rho T_ref / (gamma - 1), exact exponential over each stage.
- T_ref = c_ph^2 + psi/(n+1) inside (psi = Phi_s - Phi_wb > 0).
- Outside (superseded 10-07, see section 3): c_ph^2 for gas denser than env_amb_k (3) x the hydrostatic
  ambient density rho_amb(r,phi), c_amb^2 otherwise.
- Time constant `env_t_relax` (outside) or `env_t_relax_env` (inside), default 1e-4 code each.

### Initial state

- Discrete hydrostatic equilibrium of the WB pair, column by column. It is anchored at the deepest
  cell at the analytic n = 3 polytrope density and carried up face by face with each cell's
  isothermal background.
- This is the exact rest state of the scheme, even across the unresolved photosphere.
- Above it, the c_ph atmosphere runs until its pressure falls below the hot ambient's (or r_top);
  then the hot hydrostatic ambient at rest (section 3).
- Rotation: v_phi = (s - 1) Omega r (rotating frame) below r_top.

### Boundaries

- **r_in**:
  - `env_wall_slip = free` (default): closed free-slip wall via hydro wall_closed_ix1 (exact mirror flux:
    zero mass, energy and phi-momentum flux).
  - `noslip`: ghost v_phi held at (s - 1) Omega r.
  - Ghost density: an isothermal HSE extrapolation in Phi_wb at the edge T.
- **r_out**: as stage 1, plus energy. The stream window is at c_ph^2; elsewhere the edge T is used with
  HSE ghosts in the TRUE Roche potential (Phi_wb is flat out there).

### History `ryper.user.hst` (22 columns)

| column | meaning |
|---|---|
| Mdom, Jdom | as stage 1 (inertial Lz) |
| Menv, Jenv | mass and inertial Lz of cells with r < R_acc |
| Min, Mout | as stage 1 |
| Mwal, Jwal, Jstr | the r_in wall budget (should be 0 for a free-slip wall) |
| Jout | as stage 1 |
| MR, JR | cumulative mass and inertial Lz INTO r < R_acc through the face R_acc |
| Jin | cumulative inertial Lz in through the stream window |
| d* | the nine rates |

Budget identity: dJenv/dt = JR-rate - Jwal-rate + the tidal/indirect torque on the envelope (the phi
gravity).

### Failed attempts, kept for the record

All CPU, at quarter or half resolution, under tests2/:

1. WB everywhere (ambient included): the isothermal background of the unsupported ambient gas, with
   H_p << dr, exerts an e^7 overweight force. This gave a dt collapse at cycle 5 (tests2/dbg_a).
2. WB only on r <= R_acc (hydro/wb_rmax), with the analytic polytrope IC:
   - The top envelope cell's WB face pressure faced a plain-PLM atmosphere cell (T 728 vs 240), a
     factor ~1e3 pressure mismatch at the face.
   - This gave a NaN at cycle 223, at quarter resolution.
3. The discrete-HSE IC anchored at the TOP cell: the bottom pressure then carried column-dependent
   truncation errors, an unbalanced phi pressure gradient of ~900 (km/s)/code-time at r_in. The fix
   was to anchor at the bottom.
4. Quarter resolution (141 cells, dr = 0.05 = 15 H_p) is too coarse for the photosphere:
   - the half-cell Phi step is 7.5 c_ph^2;
   - the WB guard flattens stencils spanning more than 6 decades;
   - the result was a dt collapse at cycle 2.

   **Half resolution (dr = 0.025) is the coarsest that works**: all CPU tests below are half-resolution.

## 3. Ambient and floor redesign (10-07, after test (a) failed; user direction)

### Cause of the test (a) dt collapse (tests2/ah_env_s*, CPU half resolution, killed)

The old setup used an ambient that was cold (c_ph), uniform (rho_amb 1e-5) and at rest. That ambient
cannot support itself (H_p << dr everywhere in the domain).
- It free-falls onto the star at up to 600 km/s and piles up on the atmosphere just above r_top. Both
  are plain-gravity cells with an unresolved scale height.
- At r = 4.1234 (the first cell above r_top), etotgrav then turns the half-cell potential step (more than
  the thermal energy) into negative or garbage internal energy, and dt collapses.
- s 1.0 collapsed at t = 0.067 and s 7.2 at t = 0.066, both near phi ~ 0. The same cell went non-finite in
  stream run c_s1.0 (GPU job 12106531, t = 0.057).

### Fix

**Hot hydrostatic ambient, as the user expected.**
- Isothermal at c_amb = 300 km/s (`env_cs_amb`), at rest, hydrostatic in the TRUE Roche potential.
  Its density is `env_amb_rho` = 1e-6 at R_acc and falls by ~10x to r_out; it is floored at 1e-7.
- `rho_amb` = hydro dfloor = 1e-7 is the floor safety net: 1e-7 of the stream peak and 5e-9 of rho_ph.
- A safety-net sponge acts on gas with rho < 10 rho_amb, i.e. on the ambient itself: the velocity decays
  over `env_t_sponge` = 1e-4 and is capped at 100 km/s, and the KE removed leaves E.

**Relaxation target outside the envelope: by density relative to the hydrostatic ambient.**
- rho > 3 rho_amb,hse(r,phi) (`env_amb_k`) means stream or atmosphere material, relaxed to c_ph^2.
- Everything else is relaxed to c_amb^2.

Variants that failed (GPU, half resolution, gpu2/):

| variant | what happened | run |
|---|---|---|
| target interpolated in log rho (1e-3..1e-2) | heated cold atmosphere gas mixed into thin cells; evaporation wind, 1000-2900 km/s | a2_*, a3_* |
| two-phase attractor on T (t_mid = sqrt(c_ph^2 c_amb^2)) | heated the stream's thin leading edge to c_amb; precursor jet at 1700-1800 km/s reached the star by t = 0.0075 | g_half |
| instant relaxation (t_relax 1e-9) | the phase flip at the stream edge is a 375x pressure jump; collapse at cycle 99 | h5_inst |
| ambient density 1e-4 at R_acc | P_amb = 9 interacts with the cold atmosphere top (P up to 444 in the last cold cell): outflow at phi ~ 0 | a3_s1.0 |
| WB zone ending at R_acc (env_r_top = R_acc) | WB photosphere cell against a plain hot cell: collapse at t = 0.006 | h5b_str |

### Ambient criteria (user), measured or computed

| criterion | value |
|---|---|
| P_amb(R_acc) | 1e-6 x 9e4 = 0.09 |
| photospheric P | rho_ph c_ph^2 = 20 x 240 = 4800, ratio 1.9e-5 |
| stream ram pressure | 1 x 80.8^2 = 6.5e3 at r_out (ambient there <= 9e-3, ratio <= 1.4e-6); ~4e5 at impact |
| ambient mass, r > r_top | 3.10e-5 at t = 0, 2.85e-5 at 0.2 orbit (code mass units, h5_a1 dumps) |
| stream inflow | ~125 per code time (stage-1 hst), so ~21 per 0.2 orbit |
| ambient change vs stream | 2.6e-6 of 21, i.e. 1e-7 of the stream: negligible |
| floor-cell mass | 4.7e-6 -> 3.6e-7 |
| dt-limiting cell, envelope alone, s 1 | the ambient cell at r 4.12-4.15 (c_ad 387 km/s, \|v\| <= 60 km/s), dt 1.85e-5 at half resolution, i.e. ~9e-6 at full resolution |
| dt-limiting cell, envelope alone, s 7.2 | the deepest envelope cell at r_in (c_ad 219 + v_phi 135 km/s), dt 1.48e-5 half, i.e. ~7.4e-6 full |
| comparison | both above stage 1's stream-limited 6.58e-6 at full resolution |

## 4. Tests

All on GPU, apudev, 2 x MI300A, binary `bin/athena_ryper_gpu72_env5` (md5 1ebaa8d11b44e9cfe2e26e91d5fe518c,
commit fe16e65a source), input `bin/ry_per_accretor_env5.athinput`. Analysis scripts: scripts/ana_env.py,
scripts/dtcell.py, scripts/budget_env.py.

**Stage-1 regression.** The stage-1 input on the new CPU build vs the 1e8d0ff8 binary, 60 cycles at
quarter resolution: user.hst and hydro.hst are byte-identical (tests/bitwise_stage1.sh).

### (a) Envelope alone, no stream, 0.2 orbit, half resolution, PASS

Runs: gpu2/h5_a1 (env5, s 1, job 12108113); gpu2/a4_s1.0 and a4_s7.2 (env3 binary + env_amb_rho 1e-6,
job 12107604; s 1 matches h5_a1 to the digits shown).

**s 1.0:**
- max \|v - v0\| in r < R_acc: 1.72 km/s, in the top envelope cell (r 4.0475).
- max \|v - v0\| in r < R_acc - 0.2: 0.034 km/s.
- max \|drho/rho\|: 6.4e-4 in r < R_acc; 1.2e-5 deep.
- Menv, Jenv, Mdom, Jdom unchanged to the 6 hst digits.
- Wall fluxes are exactly 0.
- <Omega>/Omega_orb = 1.0000 at all shells (0.9999 in the top cell).
- dt: 9599 cycles, no collapse.

**s 7.2:**
- Jenv drift -2.35e-4 over 0.2 orbit. JR accounts for only 8e-11, so this is the tidal/indirect torque on
  the rotating envelope.
- Shell <Omega>/Omega_orb is 7.1985-7.1997 deep and 7.1935 in the top cell.
- max \|v - v0\|: 5.0 km/s, at the r_in wall cell (phi 5.16), with drho/rho 2.5 % there: a localized wall
  feature, not understood yet.
- dt: 11433 cycles, no collapse.

### (b) GPU smoke with the stream, full resolution 564 x 4 x 2048: smoke PASS, run FAILS at impact

**Smoke.** No FATAL in either:
- 10 cycles, job 12107709, arm b4sm (env3 + env_amb_rho 1e-6);
- env5 arm h5sm, job 12107960.

**Throughput.** 3000-4000 cycles: 2.19e8 zone-cycles/s per node, i.e. 1.10e8 per GPU.
- Measured in job 12106529 (tput) and job 12107709 (b4tp).
- That is 1.57x slower per cell than stage 1 (3.44e8): ideal EOS, energy equation, WB cache.

**Cost per orbit**, IF dt is as stage 1 (6.6e-6, i.e. 1.29e5 cycles/orbit):
- 1.29e5 x 4.62e6 / 2.19e8 = 2.7e3 s = **45 min per orbit per apu node**;
- stage 1 was 23 min.

**OPEN: the stream impact collapses dt.**
- Half resolution (env5, gpu2/h5_str, job 12107960): the stream crosses the hot ambient cleanly.
  - dt stays at 1.7e-5, limited by the ambient cell at r 4.12.
  - There is no precursor.
- dt collapses at t = 0.0499, when the stream front (rho ~1e-3, v_r -445, v_phi +465 km/s) reaches the
  photosphere at phi ~ 0.9 rad (gpu2/h5_imp, dumps every 5e-4).
  - The last dump shows dt 5.75e-6 at r 4.1234, with T = 4.5e5-5e5 in a 6e-6 to 4e-5 density cell in
    front of the impact.
  - Suspect: the WB pressure-form force on non-hydrostatic stream gas entering the WB zone (r <= r_top),
    where the isothermal background with H << dr overweights. This is the same mechanism as failed
    attempt 1, now for the stream instead of the ambient.
- Not yet fixed. Candidate fixes, not tried:
  - switch the WB force/reconstruction off per cell when \|v_r\| > some fraction of c (needs a hydro
    change, since wb_rmax is per radius only);
  - resolve the photosphere (option A cost);
  - use a hotter photosphere temperature env_cs_ph so that H_p >~ dr (a physics change).

### (c) Spin 1.0 vs 7.2 with the stream: NOT DONE

It is blocked by (b): no stream run survives the impact, so no accreted-AM comparison exists yet.

### (b) Diagnosis of the t = 0.0499 impact collapse (10-07, no code change)

Setup:
- Binary env5, input `bin/ry_per_accretor_env5d.athinput`: the env5 input plus wb_rmax = 0 (explicit) and
  an rst output at t = 0.0485.
- Base run gpu2/i_base (job 12108343): collapse at cycle 3094, t = 0.049879.
- Four restart arms from i_base/rst/ryper.00001.rst (t = 0.0485), job 12108395, dumps every 1e-5, to
  t = 0.0505. Script scripts/dtcell2.py.

| arm | change | result |
|---|---|---|
| iA | none (base restart) | collapses at cycle 3094, t = 0.049879: bitwise the same as i_base, so deterministic |
| iB | hydro/wb_rmax = 4.06: WB reconstruction and force off above R_acc (atmosphere cells 4.0726, 4.0979 plain) | **survives** to 0.0505 (cycle 3117); dt 7.2-7.8e-6 |
| iC | hydro/wb_rmax = 3.95: also the top 4 envelope cells plain | **survives** to 0.0505 (cycle 3116) |
| iD | env_t_relax 1e-5 (10x faster relaxation) | collapses EARLIER, t = 0.04960 |

### dt-limiting cell, arm iA

Up to t = 0.0496:
- It is bow-shock-heated ambient ahead of the stream front, moving in from r 4.70 to 4.12 at phi ~ 0.9.
- rho 2-16e-7, P/rho 1.7e5 -> 6.2e5 (relaxation to 9e4 cannot keep up), c_ad 530-1010 km/s.
- dt 1.2e-5 -> 7e-6.

At t = 0.04971, in cell r = 4.1234 (the first plain cell, directly above the WB atmosphere cells; r_top =
4.1106):
- rho jumps to 1.4e-5 and P/rho to 4.1e7 (c_ad 8200 km/s), with v_phi 261.
- dt 8.9e-7, then the collapse.

Arm iB has the same bow-shock ambient (P/rho 5e5) but no spike: dt bottoms at 7.2e-6.

### Confirmed mechanism

When the cold stream front reaches the photosphere, the x1 well-balanced reconstruction in the cold
atmosphere cells between R_acc and r_top is violated:
- These cells have T = 240 and H_p = 0.13 dr.
- Their isothermal background is extrapolated over a half cell by e^{+-3.75}.
- A non-hydrostatic inflow arriving there produces face states, and hence mass and energy fluxes (etotgrav
  carries the potential step), that are grossly wrong.
- These fluxes dump an enormous specific energy into the thin plain cell just above.

It is not the thermal relaxation: faster cooling makes the collapse come sooner. It is not the
bow-shock-heated ambient: arm iB survives with that same ambient.

I did not run a separate etotgrav-isolation arm. etotgrav cannot be toggled on a restart, because the
conserved energy holds rho Phi.

### The fix it points to

Restrict the WB reconstruction to the envelope, r <= R_acc (hydro/wb_rmax = R_acc). Keep Phi_wb, r_top and
the discrete-HSE IC as they are; above wb_rmax, plain PLM and the plain -rho dPhi_wb/dr force are already
coded. This needs no core change, only the existing key.

What still has to be shown:
1. Test (a) still holds with it over 0.2 orbit. The quarter-resolution, analytic-IC version of this split
   failed in attempt 2, and env_r_top = R_acc failed in h5b_str. Both differ from iB.
2. A stream run passes the impact and reaches a steady state.

If (1) fails, the per-cell alternative is a core hydro change: switch WB off where \|v_r\| > f c_s.

### Input-key fix hydro/wb_rmax = R_acc (4.06), 10-07 (coordinator go; no code change)

Binary env5, input bin/ry_per_accretor_env5d.athinput with hydro/wb_rmax=4.06 on the command line, half
resolution, GPU.

#### (1) Envelope alone, 0.2 orbit (job 12108434, gpu2/a6_s1.0, a6_s7.2): PASS, no dt collapse

| quantity | s 1.0 | s 7.2 |
|---|---|---|
| cycles | 8998 | 11436 |
| Jenv drift | 0 to 6 digits | -2.35e-4 (identical to the earlier test a) |
| max \|dv\|, deep (r < R_acc - 0.2) | 0.036 km/s (was 0.034) | 5.0 km/s at the r_in wall cell, unchanged |
| max \|dv\|, top envelope cell 4.0475 | 3.9 km/s (was 1.7) | — |
| top envelope cell <Omega>/Omega_orb | 1.069 (was 0.9999): spun up 7 % by the plain-PLM atmosphere above | 7.271 (was 7.1935) |
| atmosphere cell 4.0726 (now plain) | rearranges: dv 32 km/s, drho/rho up to 5 | dv 24 km/s, drho/rho up to 3 |
| dt-limiting cell | ambient at r 4.0979 (c_ad 387), dt 1.91e-5 | r_in wall cell (c_ad 219 + v_phi 135), dt 1.48e-5 |

In short, the deep envelope is unchanged. The top cell and the unresolved atmosphere are degraded.

#### (2) Stream run (job 12108480, gpu2/s6_s1.0): passes the first impact, collapses at t = 0.069

- The stream reaches the photosphere at t ~ 0.05 without the earlier collapse.
- dt stays 7.7-9.5e-6 (dt-limiting cell: impact-heated or stream gas at r 4.10-4.15, \|v_phi\| ~ 550-690
  km/s), then collapses at t = 0.06900 (cycle 5091; cycle 5054 in the restart).

Diagnostics:
- gpu2/s7_base: restart files at t = 0.066 (job 12108525).
- gpu2/s7_r: restart from 0.066, a dump every 1e-5 (job 12108528).

Up to the last dump before the collapse:
- At phi 1.433 (the impact longitude) the stream (rho 0.35, v_r -300, v_phi +570 km/s) runs over the
  photosphere.
- The top envelope cell (r 4.0475, WB on, rho 181, T 600) is dragged smoothly: v_phi 136 -> 141, v_r -9.
- Above it, the shock-heated plain cell 4.0726 has rho 3.6, T 3100-3200, v_phi 570, v_r -98.
- dt is steady at 7.7e-6 with no warning.

Within 1-2 cycles that top envelope cell drops from rho 181 to the floor (1e-7) with \|v\| ~ 1e21. So
this is a one-step positivity failure at the WB/plain face r = R_acc, under a Mach ~10 shear (Delta v_phi
430 km/s against c_s ~ 25-60 km/s) plus inflow. It is not gradual heating.

#### What it points to

Each of these needs a go; none was run:
- (i) hydro/fofc = true (first-order flux correction, an existing key). It needs an input-file copy,
  because the key is not in the input. I did not run it because it goes beyond the given go.
- (ii) The per-cell WB switch under inflow or strong shear (core change).
- (iii) A resolved or hotter photosphere.

(c) is still not done.

## 5. Hot-star test (user direction 10-07): c_ph = 85.3 km/s, H_p = 4 cells, WB everywhere

### Setup

Code: commit 6588011b adds env_cs_stream and env_r_hot (defaults reproduce fe16e65a).

| item | value |
|---|---|
| binary | bin/athena_ryper_gpu72_env6, md5 77a176398042496f2dafe035e7a034d1 |
| input | bin/ry_per_accretor_env6.athinput: env_cs_ph 85.3, env_cs_stream 15.5, wb_rmax 0, nghost 2, no fofc |
| grid | half resolution |
| job | 12108843 (apudev): arms t6sm, t6_a1, t6_s1 |

Choice of c_ph = 85.3 km/s:
- It gives H_p = c_ph^2/g = 4 dr at half resolution (dr = 0.0252).
- The full-resolution equivalent is 60.3 km/s (dr = 0.0126).
- The stream is injected at 15.5 km/s. Dense gas outside the star relaxes to c_ph below r_hot = r_top
  (~R_acc + 15 H_p, i.e. ~5.5 Rsun) and to 15.5 beyond.

### Physical consequence, 85.3 vs 15.5 km/s (estimates; rho_ph = 20 stream units)

| quantity | 15.5 km/s | 85.3 km/s |
|---|---|---|
| ram/photospheric pressure, 650 km/s (normal 350 km/s) | 88 (25) | 2.9 (0.84) |
| stream penetration below the photosphere | 0.027 Rsun = 8 H_p | ~0.12 Rsun = 1.2 H_p |
| envelope density contrast to r_in | 2.3e6 | 154 |
| hot hydrostatic atmosphere | — | ~14 H_p = 1.4 Rsun, T ~ 5e5 K |

**This is a numerics test, not RY Per.**

### Results

| arm | result |
|---|---|
| smoke t6sm | clean |
| envelope only t6_a1, 0.2 orbit, spin 1 | PASS, 8256 cycles, no collapse (details below) |
| stream t6_s1, spin 1, WB everywhere | SURVIVES to t = 0.31 (26608 cycles, 199 s on 1 node): the first stream run through the impact |

Envelope-only details (t6_a1):
- max \|dv\|: 0.017 km/s at r < R_acc (top cell); 0.008 km/s deep; 0.098 km/s in the atmosphere
  (rho > 1e-3 rho(R_acc)).
- max \|drho/rho\|: 8e-6.
- Jenv and Menv unchanged to 6 digits; <Omega>/Omega_orb = 1.0000 in every shell.
- dt-limiting cell: r_in (c_ad 254 km/s), dt 2.06e-5.
- That is 100-200x quieter than the c_ph = 15.5 envelope.

Stream run details (t6_s1):
- **dt history:** 2.06e-5 before the stream arrives, then 1.0-1.1e-5 after the impact (t ~ 0.05-0.07).
  The minimum is 1.00e-5 at t = 0.184. There is no collapse and no dt dip.
- **Impact point**, from the peak inward mass flux at r = 4.07 at t = 0.21-0.30:
  - phi 76.5-81.4 deg (ballistic 73.7);
  - v_r -114 to -29 km/s and v_phi 290-550 km/s (ballistic -350/+545): the stream is braked by the
    hot, dense atmosphere (rho 5-20 there).
- **Budget** (scripts/budget_env.py):

  | window | Mdot_in | Mdot through R_acc | JR / (Mdot_in j_K) | j accreted per unit mass through R_acc |
  |---|---|---|---|---|
  | t = 0.10-0.31 | 62.5 | 10.3 | 0.052 | 0.32 j_K |
  | t = 0.20-0.31 | 62.5 | 10.7 | 0.072 | 0.42 j_K |

  - Stream j_in = 1.147 j_K.
  - The instantaneous dMR and dJR change sign (an oscillating atmosphere).
  - Mdot_out is 0.4-0.7.
  - dJdom/dt is about 1.44e5, i.e. ~92 % of Jdot_in(stream) = 1.58e5. So most of the stream's mass
    and AM is piling up in the extended atmosphere or ring outside R_acc. The run is **not yet
    quasi-steady** after 0.3 orbit, and the accreted j/j_K values above are transients.
- Compare stage 1's absorbing surface: j_acc = 1.065 j_K with all of the mass accreted.

## 6. Spin pair (hot star) and option 3 (resolved real photosphere), 10-07

### Hot-star spin pair (jobs 12108898/9)

Both runs survived to t = 0.511.

**Spin 1 (t7_s1), t = 0.21-0.511:**

| quantity | value |
|---|---|
| mass through R_acc | 10.8 of the stream's 62.5 |
| AM through R_acc / (stream inflow x j_K) | 0.116 |
| j of the accreted mass | 0.67 j_K |
| j brought in by the stream | 1.147 j_K |

- In 0.05-long windows, the AM fraction through R_acc scatters from 0.015 to 0.146 and the accreted j from 0.07 to 1.19 j_K.
- About 90 % of the stream's AM is still piling up outside R_acc, so the run is not steady.

**Spin 7.2 (t7_s72): INVALID.**
- r_top (the radius above which the gas is treated as ambient) is found by searching the effective potential, which includes the spin centrifugal term. At spin 7.2 that potential peaks near corotation, ~7.5 Rsun, so the hot atmosphere never dropped by e^15 and r_top = x1max (16.24, against 6.50 at spin 1).
- The atmosphere then co-rotated at 7.2 Omega out to r_out. It was flung out (mass ~1.2e5) and choked the stream inflow from 62.5 to 13.75.
- Fix: key env_r_spin (commit ea04cd19; default r_top keeps the earlier behaviour). Staged in binary env7 (md5 a5ad0a18c1d3a73e4cb41e833f1821d5) and input env7 (env_r_spin = 4.06). The rerun is ON HOLD (user).

### Option 3: resolving the real 15.5 km/s photosphere

It is feasible with the existing polynomial + plateau radial stretch, with no core change. The user judged it too expensive, so it was not run; jobs 12109190/1 were cancelled before they started.

**Grid** (scripts/grid_design.py, arguments 900 8.0e-4 300 14):
- nx1 = 900, giving 900 x 4 x 2048 = 7.37e6 cells.
- R_acc - 0.05..+0.08 has dr = 5.9e-4..9.8e-4 Rsun (3.4-5.6 cells per H_p).
- Outside, dr/r = 2.8e-3.
- Neighbour ratio is at most 1.093.
- r_acc is set to the nearest face, 4.0596724.

**Input:** bin/ry_per_accretor_env8.athinput.

**Smoke (job 12108943):** clean.
- dt is 4.5e-7 from t = 0, set by the hot ambient's sound speed in the fine cells.
- Throughput is 2.50e8 zone-cycles/s per node.

**Cost per orbit:**

| apu nodes | cost per orbit |
|---|---|
| 1 | ~15.5 h |
| 2 | ~7.8 h |
| 4 | ~3.9 h |

These double if the impact pushes dt to ~2.5e-7.
