# FRESHGO: fresh-production smoke + launch plan (v20.3, Caltech H200)

FRESHGO, 2026-11-20 from 15:08 PST (Orion). USER target: production launch ~17:00 PST, on the coordinator's go only.

## 1. Deck changes on top of FRESHPREP (freshprep/keys + mk.sh, rebuilt)
- **hydro/eos = general, general_eos = gamma** in all decks. floor_pay = rad_column (USER, rule 1) requires eos = general (dual_energy.cpp:63). The gamma mode is the ideal gamma law through the general code path.
- **hydro/tfloor in code units** replaces tfloor_kelvin. The general EOS converts tfloor_kelvin with <units>/mu divided out (eos.cpp:318, mu_ref = 1). In gamma mode the code temperature is still p/rho with the deck mu, so tfloor_kelvin would floor at mu x T (T9: 674 K, BSG: 1959 K). The named code value is tk / (mu m_u/k_B), exactly what the ideal EOS uses:
  - BSG 423769360970.8259 (3160 K, mu 0.62);
  - m35mms 133780565496.86447 (1000 K, mu 0.6215);
  - T19/T9 (and corona) 123359972487.0939 (1000 K, mu 0.674).
- The other temperature users are mu-correct in gamma mode:
  - he_star pgen uses the <units> temperature when the table is off (he_star_m1.cpp:607);
  - rad_m1 uses temp_unit_kelvin and the tbl.active branches.
  - Only diagnostics (NanScan T print) and unused source terms (cooling, tracers) use the mu-free scale.
- **BSG inputs** read from FRESH_DIR/bsg/ (md5 = Caltech /resnick originals: 246a666e, 41edc225, e18ca805, 450bc0c1).
- **BSG blocks 256x128x128** (4 blocks; uniform full-radial x 128x128 lateral, USER/coordinator). m35mms 4 blocks, T19/T9 8 blocks (unchanged).
- **THINF variants** `*_fresh_thinf` / `*_smoke_thinf`: + rad_m1/clean_stag_thinf=true, clean_maxit=800 (addk.py). v20.3 contains v7; there is no _v key.
- Inputs: bsg-files docs/handover/fresh_prod_1120 (git 4c31e5bc0), gunzipped on Orion in freshgo/fresh_prod_1120.

## 2. CPU checks (Orion, v20.3 he build e688cb7bc md5 c5916c42, copied to freshgo/bin)
- cpu/cpu1.sh: the six smoke decks at a reduced lateral grid, 4 ranks, 2 cycles; T9 restart bitwise (4 straight vs 2 + rst + 2).
- cpu/cpu2.sh: ideal vs general/gamma equivalence (floor_pay none) for T9 and BSG, 2 cycles, bins every cycle.
- Results: see section 6.

## 3. Caltech smoke arms (bsg-files docs/handover/fresh_prod_1120/freshgo/, ARMS_smoke.txt)
- Runner armg.sh (guards on TAG/OUT) + resf.py. Every arm runs on 1 node (4 H200), NR = 4, B_v203.
- S_* arms (6 stars): a = 50 cycles straight; b = 25 + rst + 25, bitwise.
- T_* arms (THINF, 4 stars): 120 cycles, b = 60 + 60.
- B4_t9l62 / B16_bsg: block-shape timing (8 cycles); B16_bsg is also the layout-invariance check.

## 4. Production plan (fires only on the go NOTE)
Holders end **11-21 17:46 PST**, except 01-11 (testH1), which ends 11-22 00:34. Continuing past them needs the USER (OTP, new holders). Successor holders 7333072/3/35/36/37 are queued.

| star | node (holder) | layout | END (chunk wall end) |
|---|---|---|---|
| T9L6.2 | 02-11 (5840649) | 512x256x512, 8 blocks 512x128x128, 2/GPU | 2026-11-21 17:25 |
| BSG | 01-17 (hgGH3b 6133592) | 256^3, 4 blocks 256x128x128, 1/GPU | 2026-11-21 17:25 |
| T19L6.4 | 02-09 (rsgEH2 6133591) | as T9 | 2026-11-21 17:25 |
| m35mms | 01-11 (testH1 6133593) | 336x256x256, 4 blocks 336x128x128, 1/GPU | 2026-11-22 00:15 |
| T9 corona | 02-03 (hgGH3b) | as T9 | 2026-11-21 17:25 |
| T19 corona | 02-04 (hgGH3b) | as T9 | 2026-11-21 17:25 |

- Spare: 02-10 (rsgEH2), and 02-16 after RA (~17:50).
- Launcher prodg.sh STAR BIN DECK. It does a fresh start, or a restart from the newest rst (the same command resubmits).
- It sets the wall limit -t to END - 10 min, so the final rst is written before the holder ends.
- It keeps the newest KEEPRST = 3 rst files and logs chunks.txt (binary md5, deck md5, start, rc).
- Output cadence: the deck values (T9/T19 hst 1 h, bins 3 h, rst 12 h sim; BSG hst 1000 s, bins 0.5 d, rst 1 d; m35 per deck).
- Launch commands are in section 5.

## 5. Launch commands (Caltech; fire ONLY on the coordinator's go NOTE)
Common to every command:
```
G=/resnick/groups/carnegie_poc/jingze/fresh_prod_1120/freshgo     # prodg.sh + decks/
B=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_v203_e688cb7bc   # B_v203 (v20.3 e688cb7bc, md5 7df2b648); balance is off, so B_icbal would be bitwise equal
P=/resnick/groups/carnegie_poc/jingze/prod_fresh_1121              # production root (one dir per star)
F=/resnick/groups/carnegie_poc/jingze/fresh_prod_1120              # gunzipped inputs (FRESH)
MW=<caltech memwrap.sh>
```
- **Production decks (USER ~17:50 PST): `$G/decks/<star>_wngl_fresh.athinput` on B_v203 (WN-GL)**:
  - WB on, NO rebalance, GL 7x18 fixed rays;
  - THINF + clean_maxit 800 + clean_conv_max 1e-5;
  - eos general/gamma + floor_pay=rad_column, tfloor in code units;
  - 128x128 lateral blocks, paper tlim, m1_vet/m1_face at every 4th hydro_w output;
  - BSG: MESA Z0.008 tables, mu 0.62, GM 2.62502e27, 1000 K.
  - md5: bsg d4ab303f, m35mms 67d5a60f, t9l62 7935a398, t19l64 22d054ad, t9l62_corona f157df55, t19l64_corona 50105965 (bsg-files 8f922dfce).
  - WB on with the HG/RSG production keys (he_ic_balance on all but BSG);
  - clean_stag_thinf=true + clean_maxit=800;
  - fixed LQ10, 128x128 lateral;
  - eos general/gamma + floor_pay=rad_column, tfloor in code units;
  - BSG: bsg62_wb_fresh (mu 0.62, GM 2.62502e27; bsg_wb_fresh = mu 0.616 side arm only).
  - Each run stops at its deck tlim (BSG 57 d, T9 800 h, m35 10 d; T19 runs to the holder end). prodg.sh then exits, and the freed node goes back to tests (CALT3).
  - m1_vet and m1_face bins are written at every 4th hydro_w output.
- **T9L6.2:** `HOLDER=5840649 NODE=hpc-sm-02-11 OUT=$P FRESH=$F MEMWRAP=$MW END="2026-11-21 17:25" $G/prodg.sh t9l62 $B $G/decks/t9l62_wngl_fresh.athinput`
- **BSG:** `HOLDER=6133592 NODE=hpc-sm-01-17 OUT=$P FRESH=$F MEMWRAP=$MW END="2026-11-21 17:25" $G/prodg.sh bsg $B $G/decks/bsg_wngl_fresh.athinput   (mu 0.62; or the INTEG8 balance build if the USER agrees)`
- **T19L6.4:** `HOLDER=6133591 NODE=hpc-sm-02-09 OUT=$P FRESH=$F MEMWRAP=$MW END="2026-11-21 17:25" $G/prodg.sh t19l64 $B $G/decks/t19l64_wngl_fresh.athinput`
- **m35mms:** `HOLDER=6133593 NODE=hpc-sm-01-11 OUT=$P FRESH=$F MEMWRAP=$MW END="2026-11-22 00:15" $G/prodg.sh m35mms $B $G/decks/m35mms_wngl_fresh.athinput`
- **T9 corona** (test): `HOLDER=6133592 NODE=hpc-sm-02-03 OUT=$P FRESH=$F MEMWRAP=$MW END="2026-11-21 17:25" $G/prodg.sh t9l62_corona $B $G/decks/t9l62_corona_wngl_fresh.athinput`
- **T19 corona** (test): `HOLDER=6133592 NODE=hpc-sm-02-04 OUT=$P FRESH=$F MEMWRAP=$MW END="2026-11-21 17:25" $G/prodg.sh t19l64_corona $B $G/decks/t19l64_corona_wngl_fresh.athinput`
- **Resubmit** into a successor holder (7333072/3/35/36/37) after 11-21 17:46: the same command with the new HOLDER/NODE/END. prodg.sh restarts from the newest rst in $P/<star>/rst. This needs the USER (OTP, new holders).
- Every star needs a final rst before 17:30 11-21. -t = END - 10 min means the run stops by 17:15 and writes the final rst then.

## 6. Simulated time and disk budget per day of wall time
Basis: the smoke s/cycle (median, without output cycles) and dt, with the deck cadences. Rough output sizes: BSG about 3 GB per bin set in double precision (F0_bsg: 9.6 GB for 3 sets); T9/T19 about 4x the cells.

| star | dt [s] | s/cycle | sim-s per wall-s | deck tlim reached after | bin cadence -> number of outputs |
|---|---|---|---|---|---|
| BSG | 167 | 0.39 | ~430 | 57 d in **~3.2 h** | 0.5 d -> ~115 sets, ~0.4 TB |
| m35mms | 11 | 0.47 | ~24 | 10 d in **~10 h** | 3 h + 12 h -> ~80 sets |
| T9L6.2 | 97 | 1.13 | ~85 | 800 h in **~9.4 h** | 3 h -> ~270 sets, ~2-3 TB |
| T19L6.4 | 39 | 1.18 | ~33 | tlim 1e30: ~790 h of sim per day | 3 h -> ~260 sets per day, ~2-3 TB |

- **USER/coordinator decisions needed:**
  - BSG and T9 hit their deck tlim long before 11-21 17:46. Extend tlim, or let the nodes free up.
  - The total disk is several TB on /resnick (quota unknown to me). Either thin the cadences (for example, the m1_vet/m1_face bins at 4x the hydro cadence) or confirm the space.
- rst: the deck dt (T9/T19 12 h sim, BSG 1 d, m35 6 h) gives about one rst every 7-25 min of wall time; prodg.sh keeps 3. Each rst write costs ~20-75 s (smoke output cycles: BSG 22 s, T9 77 s including bins).
- These numbers will change with WB (the dt in the W_* smokes).
