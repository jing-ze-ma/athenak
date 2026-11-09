# NOTE 2026-11-09 Orion -> viper: r136fast_1109 floor-cap gates and re-block identity (answer to TASK-2026-11-09-orion-r136fast)

Short answer: gates 1, 2, 3, 4 and 6c are bitwise. 6d: L16 vs L16r1 is bitwise; L4 vs L16 is not (max rel E 6.5e-8, at the
MeshBlock seam). Gate 5: no energy is lost. The ledger residual differs between the arms by +3.4 % of dKE + dEint.
The raw Etot differs by -27 %, because the cap arm radiates that energy through the top (details below).

## Code and build
- NEW = `fastfloor-1109` 138e6ecb5 (= prod-1112-cand f34ada236 + 263ef0393 floor cap + 138e6ecb5 reblock tool).
- BASE = `prod-1112-cand` f34ada236.
- Detached worktrees `/orion/ptmp/jinma/r136fast_1109/{wt_ff,wt_base}`, kokkos submodule d8e9af03.
- Modules gcc/13 openmpi/4.1 cmake/3.30, Kokkos SERIAL. Each build took ~395 s (the two ran at the same time, -j32 each), 0 warnings:
  `cmake -B build_ff -D Athena_ENABLE_MPI=ON -D PROBLEM=he_star_m1 -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build_ff -j32`
  (same line with `build_base` in wt_base).
- Inputs: `docs/handover/consfix-1109/*` from bsg-files-1009. All MD5SUMS OK.

## Deck
`inputs/radiation/he_presn_m1_wedge.athinput` from 138e6ecb5 (identical in f34ada236). Edited with sed only, as in the consfix NOTE,
plus the TASK's mesh keys and outputs:
```
sed -e "s#HE_IC_FILE#$I/ic_he_presn_m1_gas.txt#" -e "s#HE_TABLES#$I#g" -e "s#R_TOP#2.44e11#" \
  -e '/^<mesh>/,/^<meshblock>/{s/^nx2       = 64/nx2       = 16/;s/^nx3       = 64/nx3       = 16/}' \
  -e '/^<meshblock>/,/^<time>/{s/^nx2       = 32/nx2       = 8/;s/^nx3       = 32/nx3       = 8/}' \
  -e '/^<output1>/,/^<output2>/s/^dt          = 5.0/dt          = 1.0/' \
  -e 's/^use_spherical_polar = true/&\nlayout_invariant = true\nlayout_invariant_mode = block/' $D > base.athinput
# appended: <output3> file_type=bin variable=hydro_w dt=1.0e30 bin_double=true
#           <output4> file_type=bin variable=m1      dt=1.0e30 bin_double=true
# per arm:  sed "s/^nlim       = -1/nlim       = $NLIM/"  and the cap keys after implicit_closure_thin_relax:
```
- Cap keys (in the deck, `<rad_m1>`):
  - `floorcap_rho = 9.818878e-13`. This is 30 x the minimum active rho at t = 0 (3.2729594557636787e-14, from a 0-cycle run's hydro_w bin).
    It covers 151808 of 163840 cells at t = 0 (radial index 47..639).
  - `floorcap_v = 1.0e6`, `floorcap_report = 1`.
  - floorcap_t: default (true), tmin 0 K.
  - Thermal arm: `+ floorcap_ke = thermal`.
- n_v > 0 from cycle 8 on (first `n_v 7.424e+03` at cycle 8), so floorcap_rho was not raised.
- Runs: p.shared, `--qos=serial`, one node, Slurm jobs 227903 / 227904 / 227908.
  Each arm ran as `srun -n <ranks> --exact --mem=15G --cpu-bind=cores`, all arms at the same time.
- dt is the same in the cap and no-cap arms (cycle 40 at t = 241.774 s in both, max |dt| difference over the hst 4e-7 s). On this deck, dt is not set by the floor gas.

## Comparison method
- Final hydro_w (dens, velx, vely, velz, eint) and m1 (m1_e, m1_f1..3) bins.
- Read with `vis/python/bin_convert.py` read_binary and assembled into global arrays.
- `np.array_equal` over ALL 163840 cells.

## Gates
| gate | arms | result |
|---|---|---|
| 1 keys absent, NEW vs BASE, 4 ranks, nlim 40 | G1_new vs G1_base | **PASS** bitwise (all 9 arrays, ndiff = 0) |
| 2 `floorcap_rho = 0` vs no key, NEW | G2_off vs G1_new | **PASS** bitwise |
| 3 cap on, 1 rank vs 4 ranks, nlim 40 | C1 vs C4 | **PASS** bitwise. Counters also identical: n_v 6.293504e6, n_t 1.181082e7, dKE 4.088578e32, dEint 9.728892e29 |
| 4 cap on, straight 40 vs 20 + `-r` (final rst at cycle 20) `time/nlim=40` | S40 vs C4 | **PASS** bitwise |
| 5 cap on, conservation | C4, T4 vs G1_new | **PASS** on the ledger residual: \|dResid\|/(dKE+dEint) = 0.034 (rad), 0.034 (thermal, counting the moved energy once). Raw \|dEtot\|/(dKE+dEint) = 0.27 (see below) |
| 6c re-block, nlim 0 | L4_0, L16_0, L16r1_0, S20 | **PASS**: all three are bitwise to each other and to R20's final bins. Max rel diff = 0 (all 9 arrays) |
| 6d re-block, nlim 20 | L16_20 vs L16r1_20 | **PASS** bitwise |
| 6d | L4_20 vs L16_20 | not bitwise (numbers below) |

Counter lines at cycle 39 (cumulative, all calls):
```
C4 (rad):     n_v 6.293504e+06 n_t 1.181082e+07 dKE 4.088578e+32 dEint 9.728892e+29 (erg, to E_rad)
T4 (thermal): n_v 6.293504e+06 n_t 1.181082e+07 dKE 4.088578e+32 dEint 4.098312e+32 (erg, to E_rad; KE thermalised)
S20 (c 19):   n_v 1.149696e+06 n_t 5.758976e+06 dKE 7.773454e+31 dEint 3.407634e+29
```
- In thermal mode, the clipped KE goes into e_int. The T cap then moves it on to E_rad, so dEint(T4) = dKE + dEint(C4), to 6 digits.
- The energy that reaches E_rad is 4.098e32 erg in both modes. dKE + dEint(T4) = 8.19e32 counts the KE twice.
- C4 and T4 are not bitwise to each other: dens 1.9e-9, eint 2.2e-9, E 6.7e-11 max rel. Their hst Etot agrees to about 12 digits.

### Gate 5: Etot (user hst col 19 = gas total incl. rho Phi + E_rad)
- resid = Etot(t) - Etot(0) - trapezoid of (L_in - L_top + P_esrc - P_sponge) over the hst rows (hst dt 1.0 < dt, so one row per cycle).

| cycle | t [s] | Etot no cap | Etot cap (rad) | Etot cap (thermal) | resid no cap | resid rad | resid thermal |
|---|---|---|---|---|---|---|---|
| 0 | 0.000 | 1.216485171459e39 | 1.216485171459e39 | 1.216485171459e39 | 0 | 0 | 0 |
| 10 | 71.514 | 1.227236651106e39 | 1.227236650953e39 | 1.227236650953e39 | +1.0132525176e36 | +1.0132525990e36 | +1.0132525990e36 |
| 20 | 132.769 | 1.270249316928e39 | 1.270249295644e39 | 1.270249295644e39 | +1.0159214265e36 | +1.0159231594e36 | +1.0159231593e36 |
| 30 | 188.548 | 1.336634906069e39 | 1.336634843461e39 | 1.336634843461e39 | +1.0170872600e36 | +1.0170926123e36 | +1.0170926122e36 |
| 40 | 241.774 | 1.423484712137e39 | 1.423484602575e39 | 1.423484602575e39 | +1.0194007054e36 | +1.0194148222e36 | +1.0194148220e36 |

Arm difference (cap minus no cap) against the counters (rad arm; the thermal arm's Etot and resid differences match to 9 digits):

| cycle | cum dKE+dEint | dEtot | dResid | dEtot/(dKE+dEint) | dResid/(dKE+dEint) |
|---|---|---|---|---|---|
| 12 | 1.004e31 | -2.857e30 | +2.18e29 | -0.285 | +0.022 |
| 20 | 7.808e31 | -2.128e31 | +1.73e30 | -0.273 | +0.022 |
| 28 | 1.893e32 | -5.475e31 | +4.19e30 | -0.289 | +0.022 |
| 40 | 4.098e32 | -1.096e32 | +1.41e31 | -0.267 | +0.034 |

What the numbers show:
- If the cap lost the clipped energy, dResid would track -(dKE + dEint), i.e. a ratio of -1.
- Observed: dResid is +0.02..+0.03 of the cap energy, with the opposite sign.
- The raw dEtot is -27 % of the cap energy. It equals the extra integrated L_top of the cap arm (-1.24e32 erg of the rate integral, about 30 % of the cap energy): the energy put into E_rad in the near-floor top layer escapes through the Marshak top.
- So the raw Etot difference is a trajectory difference caused by the cap (more radiation out of the top), not a loss.
- The remaining +3 % is inferred to be the hst trapezoid. Each arm's resid is ~1.02e36 (the start-up/trapezoid floor of the consfix NOTE), and dResid is 1.4e-5 of that. Energy that escapes within one ~5-6 s step is sampled poorly by the end-of-step L_top.
- I did not isolate this further (it would need a per-step flux ledger such as dbg_cons_ledger). Viper should decide whether "+0.034 ledger, -0.27 raw" meets the criterion.

### Gate 6: reblock_identity.py
```
python3 tools/e1remap/reblock_identity.py S20/rst/hepresn.00001.rst S20/hepresn.athinput rb rb20
rst t=132.76858078747838 cycle 20 nmb 4 mesh (640, 16, 16) mb (640, 8, 8) ng 3 nh 5 counts {'npred': 5, 'nt2': 13, 'neint': 1, 'nctr': 0, 'nrss': 7, 'nvgd': 0, 'ncsk': 3, 'nra': 0}
wrote rb/rb20.remap 25017120 pgen block 24
wrote deck rb/rb20.athinput        (0.9 s)
md5(rb20.remap) = 6e59260cf76f54dadc258aa751dec1a4
deck diff: + start_time = 132.76858078747838, + he_ic_remap_file = .../rb/rb20.remap
```
- R20 (6a) is the gate-4 split run S20 (NEW + cap, 4 ranks, nlim 20, final rst at cycle 20).
- Runs: `-i rb20.athinput time/nlim=0|20 meshblock/nx2=nx3=8|4` (+ `rad_m1/implicit_tol=0.9e-8` for the noise twin).
- 6d tau mask:
  - tau from the top = cumulative rho kappa_R(rho, T) dr, with T = wtemp from L4_20's rst, as in the consfix NOTE.
  - 8704 of 163840 cells have tau > 1, tau_max = 83.
  - T_wt / T_rad = 1.0000 at tau > 3.

| pair (nlim 20 from rb20) | var | all: max | all: median | tau>1: max | tau>1: median |
|---|---|---|---|---|---|
| L4 vs L16 | rho | 6.10e-9 | 2.15e-11 | 2.06e-9 | 3.62e-10 |
| L4 vs L16 | e_int | 5.65e-9 | 7.27e-11 | 2.08e-9 | 3.71e-10 |
| L4 vs L16 | E | 6.52e-8 | 1.58e-9 | 1.47e-9 | 2.30e-10 |
| noise L4 vs L4 (tol 0.9e-8) | rho | 2.55e-9 | 4.04e-13 | 6.48e-10 | 1.16e-10 |
| noise | e_int | 2.51e-9 | 1.78e-11 | 6.40e-10 | 1.26e-10 |
| noise | E | 1.33e-9 | 4.07e-10 | 6.12e-10 | 3.25e-10 |
| L16 vs L16r1 | all | 0 | 0 | 0 | 0 |

- L4 vs L16 is 2-3x the noise twin at tau > 1 (medians about equal for E).
- Over all cells, the E maximum is 49x the noise. It sits at (k, j, i) = (8, 0, 47) for both pairs. That is the L4 MeshBlock seam in x3 (k = 8) and the periodic x2 edge (j = 0), at the radius where the floor cap region starts at t = 0 (i = 47).
- At the top cells (i >= 630) the difference is 3.7e-9 against a noise of 2e-10.
- For reference only (not a gate): L4_20 (remap fresh start + 20 cycles) vs C4 (straight 40 cycles) differ by 1.3e-4 max rel in E and 8e-4 in rho, with t = 241.739 vs 241.774. The HSREMAP fresh start does not carry the dt history or the stage/predictor state, so these are not expected to agree.

## Wall times (srun step, 640 x 16 x 16, arms concurrent on one node)
| arm | ranks | cycles | wall |
|---|---|---|---|
| G1_new (no cap) / G1_base / G2_off | 4 | 40 | 263 / 265 / 263 s |
| C4 cap rad / T4 cap thermal | 4 | 40 | 194 / 195 s |
| C1 cap | 1 | 40 | 389 s |
| S20 / S40 (restart 20 -> 40) | 4 | 20 / 20 | 104 / 137 s |
| Z, L4_0, L16_0, L16r1_0 | 4/4/4/1 | 0 | 16-20 s |
| L4_20 / L4tol_20 / L16_20 / L16r1_20 | 4/4/4/1 | 20 | 136 / 137 / 172 / 329 s |

- The cap arm is 26 % faster in wall at the same dt (cpu 180 s vs 247 s on rank 0).
- BiCGStab inner iterations are about the same (527 vs 506 total), so the saving is elsewhere (not analysed).

## Files
- `/orion/ptmp/jinma/r136fast_1109/`: per-arm deck, log, hst, wall.txt; rb/rb20.athinput + tool.log; scripts (cmpbin.py, etot.py, percyc.py, dev6d.py).
- All rst, bin and rb20.remap files were deleted after the comparisons (quota).
- The worktrees wt_ff / wt_base (with builds, 0.55 GB each) are kept for reruns.
