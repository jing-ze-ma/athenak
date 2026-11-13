# TASK 11-13 viper->Caltech (USER): T9L6.2 (Jiang+2018 global LBV outburst run), everything in git (no transfer)

Files: docs/handover/t9l62-1113/ on this branch.
- Decks: t9l62_ourwall.athinput (RUN THIS) and t9l62_final.athinput (paper-style wall, optional twin).
- ICs: ic_t9l62_LIC.txt.gz and prof_t9l62_LIC.txt.gz; gunzip them.
- Opacity tables: already on this branch, docs/handover/acc3d-prod-1018/{rosseland_mesakap,planck_ext2}_gs98_x0.70_z0.02.txt (md5 d81b7a46 / 99d668f6).
- md5_unpacked.txt lists the md5s after unpacking.
- Put all four files in one dir and replace T9DIR in the deck with it.
- Full parameter table and paper sources: T9L62_SETUP_1113.md.

Binary: prod-1118 1d1972805, he_star_m1 pgen (our athena_he_gpu72_x_1d197280_p1118 / your A100 equivalent).
Run as written: 512 x 256 x 512 (r, theta pi/4..3pi/4, phi 0..pi), 8 meshblocks 512x128x128, tlim 2.88e6 s (800 h).
Cost estimate: 0.9-2.0 node-days on your nodes. Per user "follow the ic file": r_in 24.2286 Rsun, L = L_IC (implicit_flux_x1min 6.9466003350e+13).

## Inner boundary
- The user asked for OUR production inner wall (as HG/RSG): he_bc_inner=wall, he_wall_closed, he_wall_zero_flux, fixed flux, no he_bc_linhao. That is t9l62_ourwall.
- Reason: the paper-style wall (t9l62_final) injects spurious energy through (4/3) v_r E at a zero-mass-flux wall, ~2 L_wedge on day 1.
- A 1-day viper comparison of the two walls is running now. If it says ourwall is worse, an UPDATE will follow; smoke and start now.

## Do
1. Smoke 10 cycles plus a restart of 10 cycles (rc 0, bitwise restart).
2. Run to tlim. Outputs as in the deck.
3. If you have a spare allocation, also run t9l62_final as a twin.

Reply by NOTE (Caltech->viper) with smoke status, s/cycle and dt. After ~100 h, report L(r_out, t), mass-loss rate and history.
