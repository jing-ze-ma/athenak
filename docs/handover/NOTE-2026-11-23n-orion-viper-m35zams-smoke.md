# NOTE 2026-11-23n orion -> viper: M35ZAMS 256^2 S10 deck checked; Raven 8-block smoke passes exact; band memory explained

Binary: IA11 bd9553537 (A100 md5 35b66b73). Runs: Raven, 2 nodes x 4 A100 (40 GB), 200 cycles from fresh, hst only.

## (a) Deck
- We checked the deck against the M35MMS IA11 deck key by key. Only the star, grid, ray and path keys differ, plus the rst block, which is
  missing on purpose.
- All 5 IA11 keys are named and equal: time2_k1_stiff = 0, clean_stag_fwb_stage = false, clean_stag_fwb_t2 = false,
  clean_w_fse_tau = true, time2_gtime_mom = true.

## (b) Lateral band, memory, and the smokes
- **23l's "expect clamps" is wrong.** The bd9553537 default for vet_gd_raw already makes small blocks exact. The log says "raw global band
  (no clamped reads)", band 128.
- **Why the deck as is OOMs on A100.** The raw band leaves the band ring off unless vet_gd_band_ring is named (rad_m1_vetgd.cpp:1209, dcb998461
  speed-over-memory default). So the deck asks for one 42.19 GiB m1_arena per rank:
  8 B x 120 dirs x 384 shells x (128+256)(64+256).
  - The band is also not capped at the mesh extent: 114688 positions vs 65536 cells per shell.
- **The reach of 127-128 cells is a bound over all directions.** The fixed S10 set on this narrow wedge reaches only a few cells: VETREACH's CPU
  replica gives a set reach of 6 against a bound of 59. So the exact-list band of block width 64 does not clamp in practice.

| arm (8 blocks 384x128x64) | deck change | s/cycle (20->200) | dt | L_bot/L_in last [range] | L_top/L_in last | clamped reads (all ranks) | vet_gd band GB/rank (dense) | peak GPU MiB | CW |
|---|---|---|---|---|---|---|---|---|---|
| L8 | none (raw, ring off) | - | - | - | - | - | m1_arena 42.19 GiB: OOM at start | - | - |
| L8nr | vet_gd_raw = false | **1.009** | 3.4805 | 0.99989 [0.99971..1.00006] | 1.00044 | **0** | 3.21 (18.4) | 33373 | 0 |
| L8ring | vet_gd_band_ring = true | 1.328 | 3.4805 | 0.99989 [0.99971..1.00006] | 1.00044 | 0 (raw band) | 3.42 (46.0) | 29315 | 0 |

- **L8nr vs L8ring:** they are not bitwise but agree to round-off growth: L_bot 7e-8, L_in 1e-15, L_intE 1.4e-5; L_top up to 2e-3 at peak
  inside the start transient.
  - VETREACH's CPU gate finds the same: the raw band is not bitwise between 1 block and 2x2. The cycle-0 hst differs at ~1e-11 and the solver
    tolerance amplifies it. Its cause is being traced.
  - So: "exact to round-off", not "bitwise", for the raw multi-block band.
- **L4raw** (4 blocks 384x128x128, 4 A100, default keys): the exact lists hold the reach and the ring is on by default.
  - vet_gd memory: compact intensities 6.13 GB + 3 band slots 0.42 GB; ragged band 6.56 GB/rank (dense 55.2 GB).
  - It still OOMs: a 129.5 KiB allocation fails at a 39.5 GB peak, so the rest of the per-rank state for 6.3 M cells passes 40 GB.
- No dt_retry retries, stage / positivity / newton fallbacks, or FATALs in either 8-block arm.
- **Recommended keys for this fixed-ray deck on A100, until 045ed2718 lands:**
  ```
  rad_m1/vet_gd_raw = false
  rad_m1/vet_gd_allow_clamp = false
  ```
  - The second key changes the deck's current true (the [OURS-BIN] compatibility key), so any clamped read FATALs.
  - Both must be written in the deck: neither key is in it, and the command line cannot add a key.
  - This layout is the fastest and smallest: 1.009 s/cycle, 3.2 GB band.
  - vet_gd_band_ring = true is the fallback where the exact lists would clamp: exact on the raw band, 1.32x slower here.
- **Next:** VETREACH has local TEST commit 045ed2718 with the key rad_m1/vet_gd_reach = set (deck only).
  - It sizes the band from the direction set's own mu at the cell centres, uses the exact lists (no raw, no ring), and FATALs on any clamped
    read.
  - CPU gates are bitwise: set vs bound on 1 block; set (2x2, 4x4) vs the raw band; key off vs bd9553537.
  - A100 numbers (paper128 1 / 2x2 / 4x4 and this 8-block deck) will follow in a later NOTE.

## (c) Caltech H200 (141 GB)
- The 256^2 deck fits there as is, with 4 blocks of 384x128x128, one per GPU.
- Even the dense 42-51 GiB band of the no-ring raw case fits. With the default ring the band is 6.6 GB.
- The L4raw OOM says the total is above 39.5 GB. Scaling the L8 peaks (29-33 GB for 3.1 M cells) puts it at ~60 GB per GPU.
- So Caltech can run the 256^2 deck exactly now. Name vet_gd_allow_clamp = false so that any clamp FATALs.

## (d) FRESHPREP wb radii
- The rule (freshprep/keys/wb_*.txt): wb_rmax = 1.15 r_ph and he_ic_balance_rmax = 0.92 r_ph, with r_ph = the smoke-log r_int (tau = 1)
  = 5.49121e11.
- That gives **wb_rmax = 6.314892e11** and **he_ic_balance_rmax = 5.051913e11**. The deck's 6.691526e11 / 5.353221e11 do not follow the rule.
- Both are inert in this deck:
  - he_ic_balance = false, so the rmax is never read (he_star_m1.cpp:1999 / 2105);
  - wb_rmax > x1max = 5.7031e11, so the per-cell wb_rmax test never fires (hydro.hpp:830).
- Use the rule values anyway, for uniform decks.

-- orion (M35ZF)
