# NOTE from DeltaAI: faster RKG super-time-stepping (`<mhd>/rkg_lean`) and a parallel face-field halo pack

**DeltaAI, 2026-09-30 ~01:30 CDT.** Branch `sts-lean-0930`. User asked for STS resistive MHD to be closer to ideal-MHD
cost. Two changes:
- `rkg_lean`: **opt-in, default off**.
- FC pack: **always on**, bitwise-neutral.

## Measurements
All on 2 GH200. The run is 2000 cycles from the rot-302 restart of w121 scan arm `b3_e13_sts_f13` (eta cap 1e13,
RKG STS, general EOS, etotgrav, cubed sphere 24 MeshBlocks/rank). Jobs 3275840 and 3275952; scripts in the DeltaAI
tree `bench-2026-09-30-sts/`.

| path | before | after |
|---|---|---|
| STS (use_rkg_sts, s = 3) | 25.19 ms/cycle | 23.94 (FC pack) -> **21.36** (+ rkg_lean) |
| explicit Ohmic / ideal-MHD path | ~14.1 | **13.52** (FC pack) |

## 1. `<mhd>/rkg_lean` (9bf5a29b)
The RKG super-stages change only B and the total energy, and the resistive energy fluxes read only the cell-centred
field. Yet every super-stage ran all of the following:
- the full `MHD::ConToPrim`. On the cubed sphere with a general EOS that is two EOS inversions: `mhd_c2p_gen`,
  then the `cs_raisev_mhd` re-solve, ~0.8 ms per pair.
- a full conserved-variable halo exchange;
- full-array register copies.

With `rkg_lean = true`, stages `< s`:
- rebuild only `bcc0`, with `Coordinates::GnomonicCellCenteredB`. It uses bitwise the same arithmetic as
  GnomonicEquiangleRaiseVelMHD.
- skip the U exchange on both sides. `MHD::rkg_skip_u` gates InitRecv / SendU / RecvU. The Clear* functions then wait
  on already-completed (null) requests.
- copy only the IEN slot of u_ideal / u2 / uflx_ideal / u1. u1 is copied in full at stage 1.

The last stage does the full ConToPrim and then NewTimeStep. The option needs a uniform (not multilevel) cubed-sphere
mesh and is FATAL otherwise.

**Quirk of the old path, kept when rkg_lean is off:** `MHD::NewTimeStep` gates on `stage == nexp_stages` (= 2 for
rk2).
- In the RKG list it therefore fires after super-stage 2 of s, on that intermediate state.
- The `after_rkg_timeintegrator` "finalnewdt" (called with stage 1) never fires.

rkg_lean computes dt after the final stage instead.

**Checks:**
- **rkg_lean off:** the restart payload after 1 cycle is bitwise identical to the pre-change binary (2dfe4e95).
- **rkg_lean on, after 1 cycle:**
  - B is bitwise identical.
  - The total energy differs by <= 4e-13 relative. The old path repeats the etotgrav (E - rho Phi) + rho Phi round trip
    in every stage.
  - Density and momentum differ only in ghost cells: the physical BC of the last stage sees w0 one stage earlier.
  - dt differs by 4e-9, since it is computed after a different stage.
- **Over 2000 cycles:** the dt sequences of lean and base drift apart roughly exponentially (1e-5 at cycle 100, 1.5e-2 at
  cycle 2000).
- **Control:** the base run with `cfl_number` scaled by (1 + 4e-9) drifts identically (5.7e-4 vs 5.9e-4 at cycle 300,
  1.8e-2 vs 1.5e-2 at cycle 2000). The difference is the flow's sensitivity to a roundoff-sized dt change; the method
  is unchanged.

## 2. `PackAndSendFC` parallel over (m, n, v) (14047b4f)
- **Before:** the face-field pack ran 3*nmb teams (72 here), each walking all neighbours serially with a team_barrier
  after each one. It was the most expensive kernel of a cubed-sphere MHD cycle (nsys: 0.93 ms/call, 5 calls/cycle
  with STS).
- **Why the serial walk was unnecessary:** the "race condition in overlapping assignments" comment belongs to the
  *unpack*, where ghost faces are shared between buffers. The pack only reads b/cb, and each (m, n, v) writes its own
  slice `ndat*v..` of its own buffer row.
- **Now:** one team per (m, n, v), as in PackAndSendCC. The unpack is unchanged.
- **Result:** bitwise identical (restart payload vs 2dfe4e95), 0.63 ms/call.
- **Left:** it is still ~4x the CC pack (0.16 ms). Presumably the heavy seam-resample path compiled into the same
  lambda (it launches with local memory) costs the same-panel neighbours too; splitting the kernel is the next step
  if wanted.

## Not done (estimates from the nsys profile)
- **Split the FC pack:** same-panel copy kernel vs cross-panel seam kernel, ~0.4 ms/call.
- **Skip the first EOS solve in `mhd_c2p_gen` on the cubed sphere:** `cs_raisev_mhd` re-solves from the metric
  energy. ~0.4 ms/call, 4 calls/cycle; not bitwise, since the Newton warm start changes.
- These two together are estimated at ~15 % (STS) / 10-13 % (ideal) on top of the numbers above.
- **Going further for STS:** reduce per-super-stage exchanges and BCs, and fall back to explicit resistivity when
  dt/dt_diff < ~2. At e13, s = 3 is always the floor, and lean STS per sim-second only equals explicit.
