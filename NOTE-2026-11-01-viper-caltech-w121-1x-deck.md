# NOTE 2026-11-01 viper -> Caltech: WASP-121b 1x hydro fresh deck (gated) + Roche anchors (W121 1x, TOI-1518b 1x) + a p_ref issue in your 3x reference

Reply to NOTE-2026-11-01-caltech-viper-w121-r1-1x-ack (272395a5b). The files are under docs/handover/w121hfresh_1101/, next to the 3x package:
- deck/hfresh_1x.athinput (md5 800caef9...);
- README_1x.md;
- pkg/ck/CK1X/{ck,CE_tables};
- pkg/ic/ic_w121_1x_f5.txt;
- gate/, grid/;
- scripts/ana_1x.py, scripts/gate_ic_active.py;
- job_1x.sh;
- MD5SUMS_1x (written as MD5SUMS in our staging dir; paths are relative to w121hfresh_1101/, and its scripts/ana.py line means scripts/ana_1x.py).

## 1. Please check this first: p_ref in your 3x Roche reference
- Your 3x reference (ap 1.129094e10, grav 1162.669, roche_a 3.846161e11) is reproduced by the tool to 3e-6 / 5e-6, but only with **p_ref = 1.8e-4 bar**.
- That is the **10x** NRS1 transit pressure. The 3x value is **3.956e-4 bar** (w121prod_0929/met3/transit_3x.txt).
- At 3.956e-4 bar the tool gives **ap 1.141384e10, grav 1137.7664** (gate/roche_anchor_3x_pref3x.txt).
- This matters for round 2. Round 1 (1x) uses its own value below.

## 2. The 1x deck = hfresh.athinput with these keys changed (each marked `# hfresh_1x:` with the 3x value)
- `rad_met 0.0`, and also **`problem/met 0.0`, `eos_xh 0.7381`, `eos_yhe 0.2485`**: the EOS composition. These are the exact 1x-vs-3x differences of our 10-01 decks; without them the EOS stays 3x under 1x opacities.
- `ck_data_dir ck/CK1X`, `ck_table ck/CK1X/ck/Premixed_1x_g8_11_hiT2.txt`, `ck_ce_table CE_tables/FastChem_ck_1x_int_hiT2.txt`.
  - These are the same 1x tables as our 10-01 1x arms (md5 2df0d027 / b6643409).
  - **cia/, ray/ and sw_flux/ are identical to CK3X** (all 12 files checked by md5) and are not duplicated. Copy or symlink pkg/ck/CK3X/{cia,ray,sw_flux} into pkg/ck/CK1X/.
- `ic_profile ic/ic_w121_1x_f5.txt` (md5 10ee6d18...).
- **Anchor, tide off:** ap = x1min 1.122218e10, grav 1176.9609. This is planet_setup on the f5 IC with the 1x EOS, p_ref 7.07e-4 bar (1x NRS1). Our 10-01 1x values were derived on an older IC and are 0.39 % off.
- **Grid** carried over at equal pressure: x1max 2.150381e10, nx1 112, stretch keys in the deck.
  - The worst cell is 5.6 % coarser than the 3x cells/H target.
  - Band minima are within 2.3 %; at p >= 1e-6 bar the minimum is 3.09 cells/H (3x: 2.98).
- Roche off. Everything else is identical to hfresh.athinput.

## 3. Gates (viper)
- **Parse:** OK on CPU and GPU.
- **CPU fresh start, 2 cycles:** no FATAL.
  - IC T vs table (active cells): <= 3.9e-3 at 1e-9..1e-6 bar, <= 8.5e-4 at 1e-6..1e-3, <= 8e-5 deeper.
  - r(p_ref) vs the anchor: +2.3e-5 / -7.3e-6.
- **MI300A apudev, 10 min (job 12528606):** 2.57 rot, no FATAL, no eos_fail / c2p / fofc / vceil.
  - 0.0272 s/cycle; mean dt 13.05 s, 3.8 min/rot.
  - dt wanders 12.2-13.5 s instead of settling at 13.8 s as in 3x.
  - Floor events are below 3x; GTT 20.5 GiB/GPU.

## 4. Roche anchors (roche_anchor.py from c8aac483c; outputs in gate/)
- **WASP-121b 1x, tide on:** **ap = x1min 1.144050e10, grav 1132.4682, roche_a 3.846161e11.**
  - Command: `roche_anchor.py --mp 1.170 --rp 1.742 --pref 7.07e-4 --a_au 0.02571 --porb 1.27492504 --teff 6628 --rstar 1.461 --mstar 1.358 --profile ic_w121_1x_f5.txt --eos <1x EOS dump> --x1max 2.150381e10`.
  - Kepler M_* 1.394 vs literature 1.358 (+2.6 %).
  - x1max/r_L1 0.89 on the tide-off grid.
  - A grid at the Roche anchor (x1max 2.232001e10, x1max/r_L1 ~0.92, max dr/target 1.058) is in grid/gridkeys_112_roche.txt. It is **not run or gated**.
- **TOI-1518b 1x, tide on:** ap = x1min 1.257097e10, grav 1467.0475, roche_a 5.551457e11 (p_ref 1.824e-3 bar).
  - Kepler M_* +0.09 %; x1max/r_L1 0.55.
  - The Roche grid (nx1 120, x1max 2.030897e10) is not gated.
  - It goes with the TOI-1518b package in the next NOTE.

The four planet packages (WASP-76b, WASP-178b, TOI-1518b, WASP-189b) follow in one NOTE. Before that:
- WASP-76b is being re-gated at x1max 2.24e10 (0.80 r_L1, so it stays inside the lobe with the tide on).
- A 1-D UV pseudo-band test is running, because 16.5 % of WASP-178b's stellar flux (6.2 % WASP-189b, 4.3 % TOI-1518b) is below 0.26 um and band 10 deposits it ~2 dex too deep.
