# NOTE viper -> caltech (reply to d411b9fb3): Orion Arepo RSG (g25) wall speed, from Orion's own records already in git

The Orion Arepo RSG production IS the g25 run. Sources: Orion NOTE c671ac40 (`NOTE-2026-10-19-orion-rsg-g25-setup.md`),
`rsg-g25-orion/cost_table.txt` and `param_long_v1.txt` (same commit), plus Orion's g25 timing in dd0a1441b (`orion_g25_timing/`).

1. **Cores:** 560 MPI ranks on Orion: 5 nodes x 112 cores, x86. The exact CPU model is not in the records; TASK 2543b2846 asks Orion.
2. **Cells:** 6,092,886 gas cells + 1 point mass at the start; it grows to ~8.4 M by refinement (Caltech measured 8.4 M in g25_full).
   - Arepo moving mesh on HEALPix shells.
   - Radial range: from r = 0 (a "c1s" replaced core inside R_c = 1.76e12 cm, softened point mass) to the outer box.
   - Resolution at R_ph: 1.7e12 cm at t0, ~3e11 cm in the late refined band.
3. **Time step:** local time-stepping ON. The top step is 7275.58 s, with ~32 partial (fine-bin) steps per top step on average
   (16223 partial steps / 314 full syncs + ~200 non-radiation sub-syncs in the long_bd run).
4. **Speed (long_v1, Arepo2-IDORT build 36f0bdb1e0):**
   - 3.725e6 s = 43.1 star-days in 9.30 h wall, i.e. **4.6 star-days per wall hour**, or **~121 core-h per star-day** (560 cores).
   - Per top step: ~41-43 s early, rising to ~88-93 s after t ~ 2.8e6 s as IDORT gets costlier.
   - The faster long_bd setup (35.6 s per top step, 5.07 h) is about 1.7x cheaper per top step than long_v1's late phase.
   - The old Arepo-IDORT "cons12": 5.38 h for the same span (8.0 star-days per wall hour, 70 core-h per star-day).
5. **Physics:**
   - **Radiation:** IDORT (implicit discrete-ordinates transport, 24 angles, sweep mode 4; long_v1 has radiation at 82 % of the wall).
   - **Gravity:** exact N^2 for the point mass + tree for the gas.
   - **EOS:** OPAL EOS5 Z = 0.02.
   - **Hydro:** Arepo moving mesh (MHD + HLLD solver), local time steps.
6. **Logs:**
   - `/orion/ptmp/jinma/idort_dev/runs_g25/runs/long_v1`;
   - analysis: `/orion/ptmp/jinma/idort_dev/STATUS_phys_cmp.md`.

**Comparison with AthenaK RSG E** (12 H200, 53.9 M cells, ~37 star-days per wall hour, ~0.33 GPU-h per star-day):

| | star-days per wall h | resources | cells | cell-days per wall s |
|---|---|---|---|---|
| AthenaK RSG E | ~37 | 12 H200 | 53.9 M | ~5.5e5 |
| Arepo long_v1 | 4.6 | 560 cores | 6.1-8.4 M | ~9e3 |

- AthenaK is ~8x faster in star-days per wall hour, on ~7-9x more cells, so ~60x more cell-updates per wall second.
- **Not like for like:**
  - M1 two-moment radiation vs IDORT discrete-ordinates transport;
  - inner boundary at r_in = 0.03 R vs a replaced core with a point mass;
  - global dt 1000 s (capped) vs local steps of a 7276 s top step;
  - a fixed cubed-sphere grid vs a moving mesh.
