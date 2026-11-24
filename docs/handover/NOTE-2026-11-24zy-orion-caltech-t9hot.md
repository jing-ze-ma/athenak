# NOTE 2026-11-24zy orion -> Caltech: t9l62 dt: radiation signal-speed weight test (caltech-t9hot-e8eaa95af); test arms only, on a free node

Branch pushed: `caltech-t9hot-e8eaa95af` = rc3 af836f5a3 + one commit (hydro/rad_signal_speed_phase, read only when named; key off bitwise). Test arms only: no production switch (that's a USER decision after P1 / P2). Use bsg's node once bsg finishes, or any free GPUs. Never stop production for them. Note: DTREP's T is in code units (p/rho); the binders are 1e6-1e8 K, i.e. shock-hot thin gas, not a numerical heat source.


## Reading of your DTREP binders (log.5, c5150-7600, 514 lines)
- **T in DTREP is in code units, not Kelvin.** t9l62 has temp_unit_kelvin = 8.106e-9, so T_code = 1e16 is 8e7 K. The thin binders are at 8.5e5 .. 1.15e8 K, against T_rad = (E/a)^1/4 = 2.2-3.2e5 K.
  - 1e8 K is what a strong shock at the rows' own |v| gives: T_code / (3/16 v^2) has median 0.18 and never exceeds 1.
  - The gas stays hot for two reasons. kappa_P(1e7-1e8 K) is tiny: c dt kappa_P is 1e-15 .. 1e-12 in the hottest rows. And the model has no Compton term; Compton would cool this gas in ~0.4 s, against dt 10-40 s.
- **The main dt loss is the radiation signal speed, not T.**
  - In the thin rows ceff / c_gas has median 4.9, and cfl*dt_hyd / cfl*dt_gasonly has median 0.63.
  - Your ceff is reproduced to 1e-6 by rad_signal_speed_diff: c^2 = c_g^2 + (4/3 D E/rho)(1 - exp(-x)), x = tau_h c_eq/c.
  - In thin gas c_eq is 4e9 .. 2.3e10, so the radiation term grows as rho^-1/2 as the gas thins.
  - The linear radiation-acoustic root of every row moves at c_g, to < 2 %. For x << 1 the restoring part of the radiation-pressure perturbation is O(x^2); the O(x) part only damps the wave.
- The rest of the dt is kinematic. The thin rows have |v| 1-3e8 (median 1.5e8), which is 2-4x v_esc(2.5e12) = 7.7e7. Their origin is the data request below.
  - With floor_from_ekin on, the velocity ceiling (3e8; rows reach 3.02e8) keeps E fixed, so capped kinetic energy becomes heat in exactly these cells.

## New key: hydro/rad_signal_speed_phase (Orion branch t9hot-rss-rc3 e8eaa95af = rc3 af836f5a3 + 1 commit; rc4v version t9hot-rss b0ce4eef9)
- With rad_signal_speed_diff, the weight becomes 1 - exp(-x^2) instead of 1 - exp(-x).
  - At x = 1 the two weights are the same; for x > 1 the new one is larger. So dt is never larger where the radiation is trapped.
  - It changes the time step only. Without the key the code is bitwise rc3.
  - The key is read only when named. It FATALs unless rad_signal_speed and rad_signal_speed_diff are both on.
- Orion CPU results:
  - Suite: 21/21 with the key off, and 20/20 with rad_signal_speed + _diff + _phase on. That covers conservation sp/cs/lz, restart, dt_retry, layout tile, WB and order_time.
  - radwave arms: at x = 0.1, dt is 2.15x larger with c_ph error 2.7e-3 (bound 3e-3) and damping error 4.6e-4. At x = 1, dt and results are identical.
- Expected on t9l62: the thin-phase dt bound goes from median 12.5 to 22.3 (gas-only), about x1.8 on c5200-7400. It is unmeasured on the other zones.

## Arms (one test node, NOT the production node; binary = H200 CUDA 12.9 build of t9hot-rss-rc3 e8eaa95af, please push as caltech-t9hot-e8eaa95af)
Use the latest t9l62 rst of the running rc3 job (ideally one in c5200-7300). Keep all the production args exactly, plus hydro/dt_report=10. Run the four arms in parallel where nodes allow.

| arm | binary | extra args | cycles | what it checks |
|---|---|---|---|---|
| P0 | rc3 af836f5a3 (production) | - | 300 | baseline |
| P0b | e8eaa95af | - (key not named) | 50 | hst bitwise vs P0's first 50 (key off = rc3 on GPU; provenance) |
| P1 | e8eaa95af | hydro/rad_signal_speed_phase=true | 300 | dt, s/cyc, CW, FATAL, L_top, E_tot / M / L_z drift |
| P2 | e8eaa95af | the key + time/cfl_number=0.15 | to P1's end time | time-error yardstick (M35FAST hcmp.py: L_top, KE_int; P1 vs P2 vs P0) |
| P1r | e8eaa95af | the key; 100 + rst + 100 | 200 | restart bitwise vs P1's first 200 |

Post per arm:
- dt every 10 cycles;
- the DTREP binders every 50 cycles;
- s/cyc;
- the hst tail;
- the cmp results.

## Data request (from P0, or from log.5 if the outputs exist)
1. The eventlog / ecounter rows over c5100-7600 (neos_vceil, vceil_de, efloor_de, dfloor, fofc counts), or the per-cycle `<hydro> fofc_report` lines if they are on.
2. The floorcap_report lines in the same window (count, energy).
3. A sub-box around the binders, at P0's last cycle and 10 cycles earlier, as .npy. The DTREP rows of gid 0 (kji 76-79, 49-51, 58-59) and gid 5 (55-57, 122-125, 70-71) are enough.
   - Box: radial i 40-100, local j and k within +-24 of the binder.
   - Fields: rho, v1-3, T (code), pgas, E, F1-3, kappa_T, kappa_P (the opac array), the hydro_fofc count, and floorcap / vceil flags if they exist.
   - Orion uses it for the force balance of the fast thin gas: rho kappa_T (F - (4/3) E v)/c against grad p and gravity. That says whether |v| 1-3e8 is radiation driven (physical), vceil / floor-converted, or numerical.

-- orion (coordinator)
