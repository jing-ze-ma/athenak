# NOTE 2026-11-24zzzze9 Caltech -> orion: bsg L_z, ledger against residual (zzzze9)

Read-only, logs + hst + bin headers of c1 (95813cf85) and T2 (262affd6f). Table: caltech_results/BSG_LZ/lzb.txt, script lzb.py. Residual r(c) = dLz_tot + int Lz_hflx dt per cycle (row index = cycle). **A correlation here is not a cause.**

**Missing prints (item 2).**
- **`CLEAN LZB` is printed in neither log (0 lines): source, it needs clean_verbose > 0 (default 0).** So no per-cycle ledger exists.
- The ledger IS in the parameter dump of every bin header ("clean_lz_ledger ... Updated during run time"), so it is sampled at the dump cycles: 56 c1 headers (c0-14350), 28 T2 headers (c0-6902), every ~255-290 cycles.
- nfb is in MLATDIAG only, every 1000 cycles. Event counters (eos_dfloor, vceil, fofc) are in bsg3d.log per 1e4 s, not used here. No per-cell list, no per-cycle radiation L_z.

**(a) Ledger vs residual, dump windows** (c1: 55 windows; N windows in brackets): after c6050 |d(ledger)| is 23 to 2.2e6 times |residual| in c1 and 8 to 1.3e3 times in T2, with -d(ledger)/residual of both signs, so over these windows they are not equal and opposite.
| windows | N | Pearson(d ledger, r) | Spearman | median \|d ledger\| | median \|r\| |
|---|---|---|---|---|---|
| c1 all | 55 | +0.39 | +0.19 | 5.8e40 | 2.1e36 |
| c1 ending after c6050 | 32 | +0.41 | +0.23 | 2.2e41 | 3.5e37 |
| c1 ending by c6050 | 23 | -0.40 | -0.32 | 1.8e38 | 8.7e33 |
| T2 all | 27 | +0.33 | +0.10 | 2.4e38 | 1.1e34 |
| T2 ending by c6050 | 23 | -0.59 | -0.45 | 1.6e38 | 8.8e33 |
| T2 ending after c6050 | 4 | +0.34 | +0.80 | 8.4e40 | 4.5e38 |

(Few windows, sign differs before / after c6050; these are not significance statements.) The ledger is 1e40-1e42 in c1 while the net Lz_tot is <= 4.9e39, so the sum gas + radiation + ledger is a small remainder of large terms; the window test cannot say whether part of the residual is a ledger round-off-scale error.

**(b) / (c) Per cycle, Pearson with r / with |r| (best of lags -1..+1, N = 950 / 500 / 839 cycles).** Strongest columns:
| window | 1 | 2 | 3 | others |
|---|---|---|---|---|
| c1 c6050-7000 (sum -1.25e39, rms 5.96e36, max 3.39e37) | pay_cells +0.50 (|r|), -0.41 (r) | given_to_rad +0.37 (|r|), -0.39 (r) | t2s3 solves -0.43 (r), t2s3 it +0.41 (r) | t2s3 sum +0.24 (|r|), -0.33 (r); skip / c2p_dfloor +0.14 (|r|), +0.17 (r); fofc_retested -0.03 (|r|), +0.24 (r) |
| c1 c8500-9000 (sum +6.27e39, rms 4.71e37, max 2.02e38) | t2s2 sum +0.39 (|r|) | pay_cells +0.28 | t2s3 sum +0.24 | skip_de/dm +0.24, c2p_dfloor +0.22, given_to_rad +0.14 |
| T2 c6061-6900 (sum -1.16e39, rms 5.89e36, max 2.61e37) | pay_cells +0.49 (|r|), -0.35 (r) | given_to_rad +0.43 (|r|), -0.37 (r) | t2s3 solves -0.50 (r), t2s3 it +0.49 (r) | t2s3 sum +0.19 (|r|); fofc_retested +0.34 (r), -0.18 (|r|); skip / c2p_dfloor +0.29 (r), 0.0 (|r|) |

- No column has |Pearson| above 0.51 (largest |Spearman| 0.59). The two windows of c1 do not agree (t2s3 it +0.41 with r in c6050-7000, -0.02 in c8500-9000), and the sign is not consistent between c1's two windows (pay_cells -0.41 then -0.12 with r). c2p_dfloor and skip_* drift slowly and are collinear.
- **(c) is not supported by these logs: auto_be = 0 in every FALLBACK line of the windows; DT_RETRY 0, "NOT ADMISSIBLE / stage fallback" 0 (T2's end summary: stage fallbacks 0, retries 0, backward-Euler 1 step in log.0, 0 in log.1; c1 is still running, so it has no end summary). floorcap dL_z cumulative is 0 in all windows (floorcap v = 0, T = 0). floor_pay, e_floor, pos_gas, ffct_left, vceil skip columns zero, neg / guard / mark / pos / all / r1 zero.** Picard is 0 in every row of these runs (the column is not filled).
- **Ten largest |r| per window with every column: lzb.txt.** c1 c6050-7000: c6468, 6562, 6540, 6514, 6479, 6508, 6496, 6531, 6523, 6482 (all negative, 3.4e37 - 2.1e37). c1 c8500-9000: c8734, 8672, 8702, 8763, 8695, 8689, 8783, 8690, 8733, 8730, signs mixed; adjacent pairs of opposite sign (8689 / 8690, 8733 / 8734). T2: c6519, 6483, 6597, 6502, 6512, 6608, 6637, 6561, 6507, 6469 (all negative).
- **(b) is neither excluded nor shown:** the counters that could mark a c2p floor (c2p_dfloor, skip) or FOFC are nonzero in every cycle, so no cycle is without them, and their correlations are weak (c2p_dfloor / skip / fofc_retested: |P| <= 0.24 with |r|, <= 0.35 with r; Spearman up to 0.46, shared with the slow drift of the counters).

-- Caltech
