# NOTE 2026-11-01 viper -> Caltech: WASP-121b status (10-29 .. 11-01): LHLLD fix in production, deep-layer cleanup, hydro fresh start ready

FYI. Nothing to run. One build note for H100 (item 1) and the hydro fresh-start deck (item 5) are what concern you.

## 1. Production branch fork/fin51-r-fix = 182477553
- 1101f7949 (10-31): LHLLD/LHLLC odd-even fix, now the DEFAULT (item 2). WASP-121b MHD inputs need only `rsolver = lhlld`.
- 22 commits from orion (10-31): rad_m1 vet_gd cubed-sphere work and refactors (share-1031, cstwin, cslist, lb2 launch bounds, bvsplit), plus:
  - **fmad-1031**: CUDA builds with `-ffp-contract=off` now also pass `--fmad=false` (device FMA off). This affects your H100/H200 builds.
  - `-D Athena_DEVICE_FMA=ON` restores the old behaviour.
  - New binaries will not match older ones bitwise; check s/cycle once.
  - HIP (viper) is unaffected.
- 182477553 (11-01): docs/handover/scripts gets the deep-filter tools (item 3): dhj_remap.py `--hfilter-deep`, hodge_filter.py, hcv_mom.py.
- Known open: test_rad_m1_slab fails by 2.2e-8 on the 10-30 merge. It comes from the rad_m1 defaults-1030/mgsm-1030 merge, not from the solver fix. Not regenerated.

## 2. LHLLD fix (gated 10-30, merged 10-31)
- **Problem.** On WASP-121b, LHLLD left a 2-cell mode in the face-normal horizontal velocity: 1822 vs 186 cm/s for HLLD.
  - The low-Mach factor phi = O(M) removes the rho c (u_L - u_R) damping of that mode.
  - The face jump then makes the EMF numerical resistivity eta_face = (|S_M| + du/2) dx/2 negative on converging faces.
  - Result: anti-dissipation, deep equatorial upwind B term +7.2e27 erg/s vs HLLD -2.6e27.
- **Fix, both on by default:**
  - `lhlld_cb_sensor = 1` (`lhllc_cb_sensor = 1` for hydro). The 2dx fraction s of the face jump in rho u_n, from 4 cells, raises phi -> phi + (1 - phi) s. That gives HLLD/HLLC damping on a checkerboard and keeps the low-Mach benefit in smooth flow.
  - `lhlld_emf_eta_pos = true`. Clips eta_face at 0 and puts the matching heat into the energy flux.
  - Both keys at 0/false give the legacy solver, bitwise.
- **Gate (WASP g3):**
  - upwind B term -2.30e27 (sensor) / -3.05e27 (+eta_pos) vs HLLD -2.58e27;
  - 2-cell v 222 vs 186 cm/s.
  - The 2-rot pair from rot 320: the new default follows HLLD (deep eq E_mag +1.6 % vs HLLD), same cost.

## 3. Deep layers (p > ~100 bar): a fossil, now cleaned
- The deep horizontal grid-scale B is a FOSSIL from the pre-fix LHLLD run (MHD rot 300-320). The deep v checkerboard is older: hydro LHLLC from about rot 20. The fixed code stops their growth but does not remove them.
- These are NOT the cause:
  - the wall: the wall-mode test changes nothing in cells 2-6, and keeping the wall conducting is right;
  - radial resolution: halving dr in the bottom 28 cells leaves numerical resistivity unchanged;
  - horizontal resolution: C64 does not halve it in shells 2-10.
- One-time clean of the start state:
  - B: divergence-free Hodge-split filter, div B = 0 by construction.
  - rho v: shell filter with exact conservation of sum P, L_z and e + KE.
  - Arm HXV: full filter to cell 20, taper to 0 at cell 30.

| eq, 0-0.25 rot | uncleaned | HCV (filter to cell 12-20) | HXV (to 20-30) |
|---|---|---|---|
| shells 2-10 numerical/Ohmic dissipation | 2.86 | 0.36 | 0.35 |
| shells 2-10 E_mag slope [1/rot] | -0.43 | +0.014 | +0.016 |
| shells 11-20 numerical/Ohmic | 0.76 | 0.36 | 0.16-0.26 |
| shells 11-20 v 2-cell [cm/s] / odd-even | | 1058 / 0.66 | 320 / 0.18 |

- Regrowth after cleaning: +5 % per 0.25 rot in 2-10 and +21 % in 11-20. The source is the cube-panel seams (seam/interior 4-5x), a code question for long runs.
- Energy ledger and cost are unchanged by the filter.
- Production MHD start from the cleaned state is a pending decision on our side.

## 4. Smaller checks, closed
- **Deep mean-T offset** of the floor-booking fix (stswall, 5 rot): real and one-signed, but only in the 5-6 shells just below 1 bar (-0.2 K), zero deeper. It is the fix's intended cooling coming down from above, not a deep drift.
- **Mid-depth positive upwind B term under HLLD:** grid-scale induction (HLLD velocity-jump terms), not anti-dissipation.
  - The Alfvenic part takes 1.4-18x its magnetic gain from KE in shells 2-35.
  - The whole flux is dissipative globally and in shells 2-35.
  - Shells 36-50 are undecided: the per-cell booking cannot separate radial transport.
  - Budget practice: count only the field-jump part as numerical resistivity.
- **Deep convection scales** (3x IC, T_int 567 K):
  - convective below ~8 bar;
  - v_MLT 15-38 m/s, Mach 3e-3..1e-2, Rossby 0.1-0.4, eddy ~500 km;
  - resolving it needs ~C1300+ per panel;
  - the convective layer's thermal time is ~4500 yr.
  - So the global deep layer stays a quiet adiabat; a local rotating box would be cheap if needed.

## 5. Hydro WASP-121b fresh start (t = 0) on 182477553: tested on viper, GO
- Deck = the hydro twin of the s10fresh package (same grid C32 x nx1 112, ck cadence, sponges, seam keys):
  - `<hydro>` with rsolver lhllc (cb_sensor default 1);
  - **`am_vector = true`**;
  - **`problem/am_transport = false` is required**, otherwise the pgen stops (FATAL, double AM correction). The pgen then prints a harmless "L_abs drifts" warning; `<hydro>/am_transport` (default on) does the job.
- viper MI300A, 1 node:
  - 0.0277 s/cycle, dt 13.8 s, 3.7 min/rot (first 2.7 rot);
  - am_vector costs +1.5 %;
  - spin reservoir verified nonzero in every active cell;
  - 20.5 GiB per GPU; no FATAL.
- It will run on viper, so nothing is needed from you. The deck and binary recipe are in viper /viper/ptmp2/jinma/w121hfresh_1101 if you want to mirror it.
