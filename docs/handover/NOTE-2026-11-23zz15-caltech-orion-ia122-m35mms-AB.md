# NOTE 2026-11-23zz15 caltech -> orion: IA12.2 m35mms A/B (your 23zz9 steps 1-2): **B clean; s/cyc equal to A (within 0.02 %)**. Step 3 (production switch) waits for the USER's direct OK

## Build + gate
- caltech-ia122-465aff485, incremental, CUDA 12.9.41 + g++ 13.2, HOPPER90, -ffp-contract=off: **wt_orion1121/bin/athena_ia122_465aff485 md5 cfa3dc7713962f6ab2f5a1cf85de0fbe**.
- Tiny gate **PASS**: pass 1 3.686e-11, 3 passes, F/c 2.9006e+42 (caltech_results/H200CHK/gate_ia122_465aff485.txt).

## Arms
- **rst:** the newest production rst prod_ia11_1123/m35mms/rst/m35mms.00014.rst (cycle 33724, t 2.592e5, 29.1 GB), first hard-linked to rst_keep/m35mms.00014.rst, then read only.
- **Deck:** `-i ia122_decks/notile/IA122_m35mms.athinput` (md5 OK). It is **notile**, because the production m35mms runs layout_invariant_mode = block; the deck names every IA12.2 key.
- **Command line:** time/nlim=34024 (300 cycles), job/final_restart=false, output2/3/5/6/7 dt=1e30. **B** adds rad_m1/vet_gd_wall_sym=true and rad_m1/vet_gd_theta_seam_sym=true; both are reported "added from the command line".
- **Layout:** A and B ran side by side on hpc-sm-01-17, 2 + 2 H200 (CUDA_VISIBLE_DEVICES 0,1 / 2,3; 4 ranks each, 2 blocks per GPU). Max GPU memory was 73 GB.
- **Timing caveat:** the absolute s/cyc is for 2 GPUs per arm with 2 blocks per GPU, so it is not the 4-GPU production rate. Only A vs B is comparable.

| arm | CW | stage fallbacks | DT_RETRY | FATAL | s/cyc (c33750-34000) | dt first / min / last | L_top/L_in range (end) | L_bot/L_in range (end) |
|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 3.7198 | 12.0074 / 11.8881 / 11.8881 | 1.00488-1.01430 (1.00488) | 0.99188-0.99714 (0.99405) |
| B | 0 | 0 | 0 | 0 | 3.7192 | 12.0074 / 11.8550 / 11.8550 | 1.00486-1.01431 (1.00486) | 0.99188-0.99714 (0.99404) |

- A-B L_top/L_in over the 300 hst rows: rms 4.6e-5, max 1.9e-4. The dt diff key raises dt from production IA11's ~7.97 to 12.0.
- "CLEAN FALLBACK" lines (300 per arm) are the per-cycle counter prints: pos_gas 0, e_floor 0, floorcap T ~4.3e6 cells, c2p_dfloor ~100, floor_pay 0.
- **Verdict by your step-3 criterion:** B is clean (CW 0, no FATAL) and its s/cyc is within 0.02 % of A.
- **Step 3 not done:** stopping the IA11 m35mms production and continuing on B needs the USER's direct OK, which COORD is asking for.
- Files: caltech_results/IA122_M35/{A,B}/ (log_grep, user/hydro hst, log.0, run.log) and run_m35.sh.

-- caltech
