# NOTE 2026-11-21fa orion -> caltech: **GO -- fresh production (USER GO via coordinator, 2026-11-20 ~17:15 PST)**

Launch the fresh production star by star, using the commands in docs/handover/fresh_prod_1120/freshgo/FRESHGO.md section 5 (bsg-files 463a191ac).
- Decks: `$G/decks/<star>_wngl_fresh.athinput` (WN-GL).
- Binary: **B_v203 = athena_v203_e688cb7bc (md5 7df2b648)**.
- Check each deck md5 before launch: bsg d4ab303f, t9l62 7935a398, t19l64 22d054ad, m35mms 67d5a60f, t9l62_corona f157df55, t19l64_corona 50105965.
**Production has priority over every test.** Stop any test arm that holds a production node. The one exception is a WG_* smoke that still gates another star: then use any free node of the same holder (or the node that star's own WG_* just freed). **Keep R_LF on 02-16.**

## 1. BSG -- launch NOW (smoke GF_bsgnb = bsg_wngl passed: 0 CW, 0 FATAL, restart bitwise, KE 9e40)
`HOLDER=6133592 NODE=hpc-sm-01-17 OUT=$P FRESH=$F MEMWRAP=$MW END="2026-11-21 17:25" $G/prodg.sh bsg $B $G/decks/bsg_wngl_fresh.athinput`
If 01-17 is still finishing WG_m35mms, wait for it (minutes), or use another free node of a holder that ends 11-21 17:46, with END 17:25.

## 2. The other stars -- each launches **automatically the moment its WG_* smoke passes**
Pass means all of:
- 0 CW and 0 FATAL;
- restart hst + bins bitwise (rt IDENTICAL / same);
- no KE runaway, i.e. KE = hst cols 8+9+10 at the last row is below the limit in the table.

| star | smoke | limit (KE last row) | reference (WN / blown-up W) | node / holder / END (FRESHGO.md sec. 5) |
|---|---|---|---|---|
| T9L6.2 | WG_t9l62 | < 1e42 | 1.1e41 / 2.3e44 | 02-11 / 5840649 / 2026-11-21 17:25 |
| T19L6.4 | WG_t19l64 | < 3e42 | 5.5e41 / 8.6e42 | 02-09 / 6133591 / 2026-11-21 17:25 |
| m35mms | WG_m35mms | < 7e39 | W 6.6e38 | 01-11 / 6133593 / 2026-11-22 00:15 |
| T9 corona (test) | WG_t9l62_corona | < 1e42 | S 8.6e40 / W 7.8e42 | 02-03 / 6133592 / 2026-11-21 17:25 |
| T19 corona (test) | WG_t19l64_corona | < 2e42 | S 2.1e41 / W 8.4e42 | 02-04 / 6133592 / 2026-11-21 17:25 |

- If a smoke fails any criterion, **do not launch that star**; post its line and we decide.
- FRESHGO will also confirm each pass from Orion. If its confirmation and your criteria disagree, hold that star.

## 3. Report for each launch
Please give, in one commit per launch: start time, node, holder, deck md5, binary md5, and the first chunks.txt line or the first few cycle lines (dt, CW).

## 4. Holder end -- REQUIRED
- Every star must write its **final rst before 11-21 17:36 PST**. END 17:25 gives the -t stop with the final rst around 17:15.
- Then resubmit into the successor holders (7333072 / 73 / 35 / 36 / 37) with the same command and a new HOLDER / NODE / END. prodg.sh restarts from the newest rst.
- The USER may need to enter an OTP for the new holders; please flag in a NOTE when that is needed.
- Tests resume only on nodes production leaves free: v20.4c gate, then ST/STD, then H2A/H2S, then RA after R_LF.

-- orion (coordinator, for the USER)
