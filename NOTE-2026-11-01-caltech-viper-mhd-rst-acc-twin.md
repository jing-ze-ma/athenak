# NOTE 2026-11-01 Caltech -> viper (USER): MHD restart status (etafix not enough; second cause = face EMFs); cs accretor twin-off verdict

## 1. MHD restart continuation with roche_wb_tan
**On your 2a6b5e9f** (with d0d1a738):
- **The first-cycle E_B spike is FIXED:** E_B stays at 1.797e31 through t = 32 on the integ (b) STS deck. e1ecd210 still reproduces 1.8e31 -> 8.19e33.
- **The continuation pair still FAILS:** integ_b_nokeys deck, defaults on with wb_tan, 300 vs 150+150 cycles, compared from `<par_end>`. 158,223,437 of 392,956,068 payload bytes differ, the first at byte 232. You were right that it is real.

**Second cause (candidate, our probe):**
- WbTanProbe calls MHD::Fluxes on the B = 0 probe state. MHD::Fluxes also writes the six face-EMF arrays (e3x1, e2x1, e1x2, e3x2, e2x3, e1x3), which the same stage's CornerE/CT then reads.
- The probe did not save and restore them, so the stage-1 CT of the first cycle after a fresh start, and after EVERY restart, used the probe's zero-field EMF.
- **Fix:** save and restore the six arrays in ProbeT.
  - wbtan-etafix-1101 @ 267a32d0 (on top of d0d1a738, pushed; off fin51-r-fix d3f878be);
  - your integ 2a6b5e9f + 267a32d0 = 952d586b (local).
- **Being tested now on GPU:** the continuation pair (`<par_end>` compare) + the E_B check. In parallel, a second agent verifies on CPU (20+20 vs 40 cycles) and audits the probe for any other unrestored state (ck/radiation caches, RKG registers, c2p/history accumulators).

**Please also cherry-pick 267a32d0 into integ-test-1101 and keep holding the push.** We NOTE as soon as the pair is bitwise.

## 2. cs accretor: vet_gd twin OFF (USER decision), confirmed by the cs96 mass check to t 0.015

| run | max rho/rho0 (L1 tip) | cells > 20x rho0 | tip v_r |
|---|---|---|---|
| twin on + pos upwind | 2.32e4 | many (old cs96: 38622) | -3.08..+3.46 |
| twin off | 3.66 | 0 | +1.46..+2.06 (smooth, like vet_col) |

- L_meas agrees to 0.02 %. The twin removes no measurable panel/vertex imprint on cs (pole panel F_r contrasts -0.12 % both; the back-panel pattern is tidal shape).
- Twin off also looks faster (3.7 vs 5.8 s/cycle on 2 H200, not yet a clean A/B; profiling running).
- Canonical cs accretor decks: csacc-1031 @ e0f9eaef, docs/dev/accretor_rhd/cs/: twin off + twin_fuse false + twin_state mean + vet_source_noesrc false, implicit_pos_transport upwind, and problem/floor_m1_vcap true in the stream decks (USER-approved).
- **vet_gd_twin_srw** (the bitwise-safe fix for runs that keep the twin) is in csacc-1031 as an option; nothing uses it.
- **No global twin default was changed.** HG/RSG keep the twin (thick-cell check clean). An RSG twin-off A/B (25 sim days, one node) is running to see whether production could drop it too.
