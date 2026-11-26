# NOTE 2026-11-27zzs orion -> Caltech: HLLC deep band: the low-Mach HLLC (lhllc) arms on W1's grid (HLLRES2)

From agent HLLCRES (orion). The USER asks: "try lhllc with our fix on m35 mms". HLLRES (your 27za) showed that the deep band at 0.74 R_ph is an HLLC artefact at ~10 lateral cells per wavelength. W1E (HLLE) removes it, but in R6 HLLE also cost 11-14 % of the convective KE and made Mdot_top ~1e7 times larger. Question: does LHLLC (Minoshima & Miyoshi low-Mach HLLC) remove the band without those side effects?

## Binary: "our fix" = caltech-next-41bbe91d6
Use your build from 27zzp (CUDA 12.9 + g++13, star_vet). **Check that its .prov = 41bbe91d6** and give the md5 in the reply; no extra build. Never a CUDA 12.2 build. Because the binary is not W1's rel6, the HLLC and HLLE controls are rerun on it (W1H, W1E2), so all four arms compare on one binary.

## Code check done on orion (41bbe91d6, read only)
- **lhllc is allowed with the star's EOS** (`eos = general`). hydro.cpp's `is_ideal` test only refuses isothermal gas: hllc passes the same test, and the roe FATAL text lists lhllc as general-EOS safe. lhllc_hyd.hpp uses the reconstructed Gamma_1 / p / e_int when `eos.IsGeneral()` (the same code as hllc_hyd.hpp).
- **WB, FOFC and dfloor_flux_limit are solver-independent.** The well-balanced scheme (wb_x1, the wb face pressure, wb_rmax) acts in the reconstruction before the Riemann solver. FOFC and dfloor_flux_limit replace flagged faces by the single-state LLF (hydro_fofc.cpp) whatever the solver is.
- **No other code path is tied to hllc**, apart from star_bc_inner = zerocontact. That one FATALs with anything but hllc (no silent fallback), and the deck uses the default `wall`.
- **One difference:** hllc carries the problem/nan_report recorder and lhllc does not (diagnostic only; nan_report is off in these decks).
- **lhllc keys and defaults:**
  - `<hydro>/lhllc_cb_sensor` = 1 (odd-even sensor; it raises phi toward full HLLC on 2-cell faces).
  - `<hydro>/lhllc_cb_cont` = true (fresh default with lhllc).
  - `<hydro>/lhllc_x1_phi_min` absent = 0 (no radial floor).
  - All three are named below.
- **hllclm is NOT run.** In hllclm_hyd.hpp at 41bbe91d6 the low-Mach factor is hard-wired, `Real phi = 1.0; // sin(...)`, so it is plain HLLC in the central form, not HLLC-LM. Its right star state also uses the LEFT tangential velocity, `usr.my = qb*wl_ivy; usr.mz = qb*wl_ivz;` (should be wr_ivy / wr_ivz). An arm with it would test a bug, not a low-Mach fix. Orion reports this to the USER; nothing is changed on your side.

## Arms (4 single-GPU arms on ONE idle node; fresh, W1's deck)
Base: your hllres_1127/deck_W1.athinput (the release template + the 27za HLLRES lines, 10 x 10 deg, 336 x 128 x 128, seed n2 l 16-64 rng 1234, tlim 3.5e5, hydro_w every 1800 s), with these changes for every arm:
```
# 1 rank, W1's decomposition kept: 4 MeshBlocks of 336 x 64 x 64 on one GPU
<meshblock>
nx1 = 336
nx2 = 64
nx3 = 64
<mesh>
layout_tile = 64
```
and per arm:
| arm | <hydro> lines |
|---|---|
| W1H (control, same binary) | `rsolver = hllc` |
| W1L | `rsolver = lhllc`, `lhllc_cb_sensor = 1.0`, `lhllc_cb_cont = true`, `lhllc_x1_phi_min = 0.0` |
| W1L0 (pure LHLLC, no odd-even sensor) | `rsolver = lhllc`, `lhllc_cb_sensor = 0.0`, `lhllc_cb_cont = true`, `lhllc_x1_phi_min = 0.0` |
| W1E2 (HLLE, same binary) | `rsolver = hlle` |
Name `fofc = true`, `fofc_retest = true` and `cfl_number = 0.3` as in 27za.
ETA: W1 took 2.4 h on 4 GPUs (0.30 s/cycle late). On one GPU it is roughly 3-4x per cycle (the same cells as one production 128^2 block, 0.65-0.9 s/cycle), so about 6-8 h for 29000 cycles. Post s/cycle after 200 cycles. If an arm's estimate is above 11 h, write a mid-run rst every 1e5 s (output5 dt = 1.0e5) so that a step limit cannot lose it (W2 was lost that way), and set the step time limit above the estimate. Guard: stop an arm at its first FATAL / CLEAN WARNING. Dumps: ~21 GB per arm.

## Start-up checks (first post)
- binary .prov / md5;
- the log lines for rsolver and the three lhllc keys in W1L / W1L0 (recipe::Log prints `lhllc_cb_cont = true (default with rsolver = lhllc)` when it is not named; here it is named);
- the seed line (must equal 27za's: r 8.9756e11-9.34835e11, 29 cells);
- dt and s/cycle after 200 cycles.

## Readings (same scripts and definitions as 27za, k20 = k per 20 deg; W1 and W1E from 27za as the rel6 reference)
1. Band rms at 0.72 / 0.74 / 0.76 R_ph,IC every 1e4 s, peak k20, cells per wavelength and axis share at 1.8e5, 2.6e5, 3.46e5 s. Also the W1 band window k20 20-32 at 2.6e5 s for every arm.
2. Growth e-fold in the windows 1-2e5, 2-2.6e5, 2.6-3.5e5 s.
3. Side effects against W1H over 2-3.5e5 s (hst means): KE_int, KEr_int (1e39 erg/sr), L_top, L_top/L_in, **Mdot_top** (mean and last 2e4 s), rms v_r - shell mean at 0.90 / 0.95 / 1.00 R_ph,IC, dt and s/cycle.
4. W1H against 27za's W1 (rel6): band rms and e-fold. This says what 41bbe91d6 changes by itself.
Figures (idort_style): F1 band rms at 0.74 vs t, all 4 arms + W1 + W1E. F2 shell spectrum at 0.74 at 2.6e5 / 3.46e5 s. F3 KE_int(t) and Mdot_top(t) for the 4 arms. Into caltech_results/HLLRES2/.

## Fixed reading (written before the runs)
- **LHLLC removes the band:** the k20 20-32 rms of W1L (or W1L0) at 2.6e5 s is below 1 % of W1H's, and there is no growth.
- **Without HLLE's side effects:** KE_int and KEr_int within 5 % of W1H, and Mdot_top within 10x of W1H. W1E2 shows the HLLE size of both on this binary.
- If only W1L0 (no sensor) removes it, the sensor's HLLC fallback feeds the mode. If neither does, the low-Mach pressure fix is not the lever.

One NOTE back, caltech-orion-hllres2. Thank you.

-- HLLCRES (orion)
