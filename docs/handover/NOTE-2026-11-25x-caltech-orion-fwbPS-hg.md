# NOTE 2026-11-25x caltech -> orion: prototype 4a7885d73 on HG (cubed sphere), 100 cycles, T6 (your 25x)

## Result (FW = 4a7885d73, NX = 79c1cabef; no dcycle dumps, only the driver's final bins)
- FW: 2026-11-25 04:44:34 rc=0 ended: Terminating on cycle limit wall=327 s CW 0 FATAL 0
- NX: 2026-11-25 04:50:06 rc=0 ended: Terminating on cycle limit wall=332 s CW 0 FATAL 0
| arm | rows | s/cyc mean / median (elapsed, first and last dropped) | dt min/mean/max | its per cycle | CW/FATAL/inf/DT_RETRY | fallback / BE |
|---|---|---|---|---|---|---|
| FW | 100 | 2.8090 / 2.7613 | 90.000/90.000/90.000 | 30.8 | 0/0/0/0 | 0.000000e+00 / 0.000000e+00 |
| NX | 100 | 2.8226 / 2.7735 | 90.000/90.000/90.000 | 30.6 | 0/0/0/0 | 0.000000e+00 / 0.000000e+00 |

| cycle | t | dt FW/NX - 1 | L_top rel | L_bot rel | L_mid rel | L_int rel | L_midE rel | L_intE rel | KE_int rel | Lz_tot rel | Lz_hflx rel | M_tot rel | Etot rel | Mdot_top rel |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 6.609609e+07 | 0.00e+00 | 4.711e-12 | 8.231e-10 | 5.125e-02 | 4.780e-02 | -7.292e-06 | -2.096e-05 | 1.370e-10 | 2.220e-16 | 8.825e-02 | -1.110e-16 | 2.220e-16 | 1.928e-02 |
| 10 | 6.609690e+07 | 0.00e+00 | 1.091e-05 | -1.223e-09 | 4.496e-02 | 3.743e-02 | 1.827e-05 | 7.156e-05 | 3.497e-09 | -2.109e-15 | -1.184e-02 | -8.882e-16 | -5.214e-12 | 5.438e-03 |
| 50 | 6.610050e+07 | 0.00e+00 | 1.090e-04 | 4.188e-10 | 2.859e-02 | 2.784e-02 | -2.006e-04 | -3.792e-04 | 1.893e-08 | 9.326e-15 | 9.670e-03 | -6.661e-16 | -1.697e-09 | 6.380e-03 |
| 100 | 6.610500e+07 | 0.00e+00 | 1.853e-04 | 1.141e-11 | 3.838e-02 | 3.570e-02 | 2.487e-04 | -3.854e-04 | 3.704e-08 | 1.199e-14 | -4.084e-01 | 3.109e-15 | -4.903e-09 | -1.071e-01 |
- L_top / L_in: FW c1 0.963075, c10 0.971615, c50 0.976093, c100 0.973498; NX c1 0.963075, c10 0.971604, c50 0.975987, c100 0.973318
- values at the last row, FW | NX: L_top 1.36646127e+38 | 1.36620818e+38; L_bot 1.50545081e+38 | 1.50545081e+38; L_mid 1.36026832e+38 | 1.30999400e+38; L_int 1.36174187e+38 | 1.31479847e+38; L_midE 1.35931754e+38 | 1.35897958e+38; L_intE 1.35797766e+38 | 1.35850121e+38; KE_int 1.92851214e+44 | 1.92851207e+44; Lz_tot -7.77901999e+47 | -7.77901999e+47; Lz_hflx 5.23325463e+31 | 8.84532151e+31; M_tot 7.07138780e+31 | 7.07138780e+31; Etot 3.63193105e+46 | 3.63193106e+46; Mdot_top 1.46638590e+15 | 1.64223677e+15
- FW CLEAN FCONE lines 200; first: CLEAN FCONE cycle 794784 t2s 2 faces 428119 momentum cut sum |dF| V/c^2 1.280117e+27 | last: CLEAN FCONE cycle 794883 t2s 3 faces 803731 momentum cut sum |dF| V/c^2 2.881919e+27
- FW CLEAN FWB lines 200; first: CLEAN FWB cycle 794784 t2s 2 faces 150026140 radiation momentum sum |dp| 4.955995e+29 x1 faces sum dp 2.482113e+28 | last: CLEAN FWB cycle 794883 t2s 3 faces 150095994 radiation momentum sum |dp| 6.643671e+29 x1 faces sum dp 1.404205e+28
- NX CLEAN FCONE lines 200; first: CLEAN FCONE cycle 794784 t2s 2 faces 36235796 momentum cut sum |dF| V/c^2 1.557518e+29 | last: CLEAN FCONE cycle 794883 t2s 3 faces 29938565 momentum cut sum |dF| V/c^2 4.706697e+28
- NX CLEAN FWB lines 200; first: CLEAN FWB cycle 794784 t2s 2 faces 151577909 radiation momentum sum |dp| 5.092344e+29 x1 faces sum dp 8.701601e+27 | last: CLEAN FWB cycle 794883 t2s 3 faces 151645029 radiation momentum sum |dp| 7.350602e+29 x1 faces sum dp 1.622151e+28
