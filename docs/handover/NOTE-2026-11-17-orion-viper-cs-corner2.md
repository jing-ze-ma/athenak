# NOTE 2026-11-17 Orion -> viper: cs corner follow-ups (USER-approved): stratified hydro + floor, and M1 corner

Follows NOTE-2026-11-17-orion-viper-cs-corner.md (hydro blob test: no key acts, no corner jet).

## Answer
- **No key points at the RSG E corner jet in these CPU tests either. No GPU twins recommended.**
- **Stratified atmosphere with floor (1e10 contrast), gravity, WB and merge_lat:**
  - guard 4 vs 0 is **bitwise identical**, with and without merge_lat.
  - `cs_vertex_3cell` acts there (FOFC fires in the floor), but changes the corner metrics by < 1 %.
  - Corner, seam and panel-centre max |v_r| agree within 1-3 % in every band, up to t = 6 (22 sound-crossing times of H).
- **Hot buoyant bubble rising into the floor layer at a vertex vs a panel centre:**
  - guard 4 vs 0 is again BITWISE.
  - The corner plume is about 2 % weaker on axis and about 10 % wider laterally. That is not a jet.
- **M1, optically thick hot spot:**
  - The keys act, near the vertex. Production vs all-off (`cs_vertex_mono`, `cs_cross_mono`, `cs_vertex_3cell`) moves the corner state by at most 2 %.
  - Corner/centre ratios of E, F_r, T, rho, v_r are 1.00-1.03 in both arms. No focusing, no off-axis jet.
- **Both restart checks (USER rule) are BITWISE** (data arrays; headers differ only in the echoed parameters):
  - stratified deck (WB + merge_lat + FOFC + guard), 40 straight vs 20 + rst + 20;
  - M1 production deck, blob at the corner, 16 straight vs 8 + rst + 8.

## Code / binaries
- Fork `race3-prod-1118-viper` 0e7f5ca26 (CPU, `-DPROBLEM=cs_test`, MPI), plus a local throwaway pgen patch: /orion/ptmp/jinma/cscorner_1117/blob_strat.patch, **not pushed**. The patch:
  - adds the 3-D blob option to iprob 12 and 13;
  - for iprob 13 fills `phicc0` / `phi0` with phi = g (r - r0) so `<hydro>/wellbalance_dynamic` sees the potential (constant g, applied by CSTestGravSrc as before);
  - re-enrols the iprob 13 BC/source/potential on a restart. Before this patch, cs_test returned early on restart, and a restart of iprob 13 died with "User BCs ... not enrolled".
- Stratified binary: /orion/ptmp/jinma/cscorner_1117/wt/build_bub/src/athena (restart check, bubble); build_strat (same source minus the restart and bubble edits; main arms).
- M1 binary: /orion/ptmp/jinma/cscorner_1117/wt/build/src/athena (no change; `rad_m1` is created when the `<rad_m1>` block exists).

## 1. Stratified cs atmosphere with floor (decks /orion/ptmp/jinma/cscorner_1117/strat/base*.athinput)
- **Atmosphere:** cs_test iprob 13, isothermal and hydrostatic.
  - Constant g; r from 1 to 2; H = 0.0434, so rho falls by **1e10** over the domain.
  - dfloor 1e-9: the top ~2 H sit at the floor. pfloor 1e-13, dfloor_keep_velocity.
  - Cell noise vnoise 1e-3 c_s. User radial BCs hold the exact state.
- **Grid:** N = 32 per panel, nx1 = 96 (4.2 cells per H), 2x2 MeshBlocks per panel.
- **Hydro keys as in the RSG deck:**
  - plm, hllc, fofc, nghost 3;
  - `wellbalance_dynamic`, `wb_x1`, `wb_option = polytropic`;
  - `cs_seam_flux = shared`, linear resample, blockend extend, `cs_vertex_3cell = true`;
  - am_transport on (default).
- **merge_lat arms:** `merge_lat_rc = 1.1:1.2:1.35`, `merge_lat_f = 8:4:2`, rigid, `blockend = exchange`.
- General EOS not used (needs tables; not trivial in cs_test).
- tlim 6 (about 22 H/c_s, 0.7 lateral panel crossings), 8 MPI ranks each on hold 230761.
- **Arms:** guard 4 / 0, each with and without merge_lat. Plus `cs_vertex_3cell = false` (with guard 4 + merge_lat, and with guard 0 without merge_lat); the 3cell arms were stopped at t = 3.6 / 2.4 for footprint and compared at t = 1.5 and 3.0.
- **Do the keys act?**

  | comparison (t = 6.0 unless noted) | result |
  |---|---|
  | guard 4 vs 0, merge_lat | BITWISE (all of rho, v, e) |
  | guard 4 vs 0, no merge_lat | BITWISE |
  | merge_lat on vs off (guard 4) | max d rho 1.6e-4, d v_r 3.4e-2 (the merged shells change the flow) |
  | 3cell off vs on, t = 1.5 / 3.0 | not bitwise: d v_r 8e-3 / 3.9e-2, all in the floor layer (rho0 ~ 1e-9, r ~ 1.9). Above the floor (rho0 > 1e-7) d v_r is 1.8e-6 / 2.4e-4. The differences start ~0.15 rad from a vertex and then spread over the sphere |

  - FOFC fires (fofc_report): ~1e7 cell-stages per 1.5 time units, 98 % of them at the density floor. That is why 3cell acts here and not in the unstratified test.
- **Corner vs seam vs panel centre:** max over cells within 0.08 rad of the 8 vertices / 12 seam midpoints / 6 panel centres, guard 4 + merge_lat (the no-merge_lat arm is the same within a few %).
  - Bands by initial density: deep rho0 > 1e-3, mid 1e-7 .. 1e-3, top < 1e-7 incl. floor.

  | t | band | corner max abs v_r | seam | centre | corner max d rho/rho | seam | centre |
  |---|---|---|---|---|---|---|---|
  | 1.5 | deep | 3.03e-3 | 3.01e-3 | 3.00e-3 | 2.8e-4 | 2.5e-4 | 2.4e-4 |
  | 1.5 | mid | 3.83e-3 | 3.80e-3 | 3.79e-3 | 2.7e-3 | 4.2e-3 | 3.5e-3 |
  | 1.5 | top | 0.291 | 0.291 | 0.290 | 0.143 | 0.139 | 0.112 |
  | 3.0 | mid | 8.10e-2 | 8.02e-2 | 7.90e-2 | 4.2e-2 | 3.9e-2 | 4.4e-2 |
  | 3.0 | top | 0.263 | 0.262 | 0.260 | 0.76 | 0.89 | 0.63 |
  | 6.0 | mid | 9.27e-2 | 9.05e-2 | 9.02e-2 | 1.9e-2 | 2.2e-2 | 2.2e-2 |
  | 6.0 | top | 0.223 | 0.222 | 0.216 | 0.093 | 0.092 | 0.065 |

  - The floor layer runs at ~c_s everywhere (top band max abs v_r ~ 0.22-0.35 with c_s = 0.27). This is the floored, non-hydrostatic top, and it is the same at corner, seam and centre.
  - **Density imprint:** in the top/floor band the density deviation at the vertices and seams is 1.2-1.4x the panel-centre value. That is a seam/corner imprint in the floor, but it carries no velocity excess.
  - **v_r:** corner/centre is 1.01-1.03 in all bands and times.
  - With merge_lat, the location of the global max abs v_r below the top band is 0.06-0.10 rad from a vertex at t = 3-6 (r ~ 1.69, just below the top band), but it exceeds the centre-region max by only 2-3 %. Without merge_lat it is 0.25-0.30 rad from a vertex.
- **Restart check** (rstchk/): base deck with merge_lat + guard 4, 8 ranks. 40 cycles straight vs 20 + restart + 20: **BITWISE** in every hydro_w array at cycle 40. rst files deleted.
- **Hot buoyant bubble in the same atmosphere** (lateral density contrast across seams, the regime the guard was built for):
  - blob p x31 (T x31, top hat, w 0.15) at r = 1.35 (rho0 = 3e-4), placed at C (1,1,-1) or P (0,0,-1), guard 4 vs 0, t = 0.5. The seam arms were stopped for footprint.
  - **guard 4 vs 0: BITWISE** at both placements, so the guard never fired. The lateral contrast stays below 4x across the seam resample.
  - The plume reaches the floor layer (r ~ 2) in both placements by t = 0.5. Both are shown as C / P:

    | quantity | C | P | C/P |
    |---|---|---|---|
    | axis peak v_r | 1.41 | 1.44 | 0.98 (corner weaker) |
    | max v_r below the floor band | 0.90 | 0.87 | 1.03 |
    | v_r in the 0.15-0.30 rad ring around the axis (r 1.3-1.6), max | 0.75 | 0.67 | 1.12 |
    | same ring, median | 0.39 | 0.35 | 1.11 |

  - At the corner the plume is about 10 % wider laterally and about 2 % weaker on axis: a mild defocusing, not a jet.
  - No non-finite cells.
- Plot: docs/handover/cscorner-1117/cs_corner_strat.png (max abs v_r and d rho/rho vs t, corner/seam/centre, all four arms).

## 2. M1 radiation corner test (a worker on Orion; deck /orion/ptmp/jinma/cscorner_1117/m1/base.athinput)
- **Grid and hydro:** cs_test iprob 12 hot blob (p and T x4, w 0.2, r 1.6), N = 48 per panel, nx1 = 48 over r 1-3, 2x2 MeshBlocks per panel. Hydro as above.
- **Radiation:** `<rad_m1>` = the production block of docs/dev/accretor_rhd/cs/canon_cs96_static_s1.athinput (implicit transport, vet_col + vet_gd, hesdirk2, blend idort_f, ...), with:
  - constant opacity kappa = 30 (tau ~ 1 per cell, optically thick);
  - c = 30 (~20 c_s), arad 0.1;
  - vet_scatter off, force_reference none;
  - reflect / marshak radial BCs.
- tlim 0.2 (~41 cycles). E_rad starts at 0 and reaches aT^4 in the first step, so the blob becomes a radiation-dominated hot spot.
- **Arms:** production (`cs_vertex_mono = true`, `cs_cross_mono = true`, `cs_vertex_3cell = true`) vs all-off. Positions C (1,1,-1), E (1,0,-1), P (0,0,-1).
  - The single-key arms were started at the corner only, with tlim 0.1. They were stopped for footprint (USER: corner follow-ups are lower priority), so there is no single-key attribution. The one that finished (cross_mono off) has its dump at a different time from production, so it is not comparable.
- **Do the keys act?** Max abs(all-off - production) relative to max abs(field):

  | pos | t | rho | v_r | E | F_r | F_lateral | where |
  |---|---|---|---|---|---|---|---|
  | C | 0.1 | 3.2e-3 | 5.3e-3 | 3.8e-3 | 2.1e-2 | 6e-2 | within 1-3 deg of the (1,1,-1) vertex; ~100x larger there than > 30 deg away |
  | E | 0.2 | 6.9e-4 | 8.8e-4 | 1.6e-3 | 3.5e-3 | 1-5e-2 | at the cube vertices |
  | P | 0.2 | 6.5e-5 | 7.6e-5 | 9.3e-5 | 3.0e-4 | 0.7-5e-3 | far from the blob |

- **Anisotropy:** peak along the blob axis of the excess over the same-radius shell median.

  | t | arm | pos | E | F_r | T | rho | v_r |
  |---|---|---|---|---|---|---|---|
  | 0.2 | prod | C | 0.0667 | 0.0534 | 0.1586 | 0.0653 | 0.0867 |
  | 0.2 | prod | E | 0.0673 | 0.0569 | 0.1597 | 0.0677 | 0.0909 |
  | 0.2 | prod | P | 0.0663 | 0.0526 | 0.1578 | 0.0658 | 0.0862 |
  | 0.1 | prod | C | 0.153 | 0.176 | 0.292 | 0.0313 | 0.170 |
  | 0.1 | all-off | C | 0.153 | 0.174 | 0.292 | 0.0308 | 0.169 |

  - Corner/centre at t = 0.2 is 1.006 / 1.014 / 1.005 / 0.992 / 1.006 (E, F_r, T, rho, v_r).
  - At t = 0.1, corner/centre is 1.01-1.03 in production and 1.006-1.02 in all-off. The production keys make the corner excess slightly *larger* (+1.6 % rho, +1 % F_r) than all-off.
  - Seam/centre is 1.01-1.08, slightly larger than corner/centre.
  - The fronts (10 % of the axis peak) are at the same radius for all positions.
  - The highest v_r at r > 1.9 is 0.1377-0.1379 in every run, on the blob axis. There is no off-axis or corner jet.
  - With the blob elsewhere, the excess in the corner direction is 1e-4 (P) to 4e-3 (E), so the vertex creates nothing on its own.
- **Restart check:** production, blob at the corner, 16 cycles straight vs 8 + rst + 8. All hydro_w and m1 E, F1-3 arrays are BITWISE. rst files deleted.
- Plot: docs/handover/cscorner-1117/cs_corner_m1_axis.png (the axis profiles of C/E/P production and E/P all-off overlap).

## What is still not covered
- General EOS and tables, real (non-constant) opacity, the optically thin regime.
- Convection: these tests have no convective driving.
- Long times. The RSG jet built up over hundreds of days; these tests run for tens of dynamical times.
- **Single-key attribution in M1.** If wanted, rerun production and each single-key-off arm with the same nlim (about 20 cycles), so the final dumps are at the same time.
- The CPU tests cannot reproduce the production event itself. They show only that none of these keys makes the cube vertex focus or launch flow in isolated setups:
  - hydro;
  - stratified floor with WB and merge_lat;
  - thick M1.
