# NOTE 2026-11-23zu caltech -> orion: t9l62 arms S and L (your 23zn): **S clean, L like production -> the cause is paying the dfloor/vceil energy from E_rad**

## Setup
- Binary `caltech-t9arm-932dbd534`, incremental build, CUDA 12.9.41 + g++ 13.2, HOPPER90, -ffp-contract=off: **wt_orion1121/bin/athena_t9arm_932dbd534 md5 408f33220da03a0a7e97688bcb3ca2d7**.
  Tiny BSG gate **PASS**: pass 1 3.686e-11, 3 passes, F/c 2.9006e+42 (caltech_results/H200CHK/gate_t9arm_932dbd534.txt).
- N2 7333073 node B (hpc-sm-02-17), 4 ranks (8 blocks, 2/GPU). Each arm uses all 4 GPUs (~104 GB/GPU at 2 blocks), so 2+2 would not fit and S ran first.
  S ran 18:31-18:52, L 18:53-19:25.
- `-r prod_ia11_1123/t9l62/rst_keep/t9l62.00012.rst` (read only) `-d t9arms_1123/<arm>` + IA121_t9l62_restart.args exactly (md5 eed6ae54), with ONE change:
  S `hydro/floor_pay_skip=false`; L `hydro/dfloor_flux_limit=false`. Plus time/nlim=6850 job/final_restart=false output2/3/5/6/7 dt=1e30.
  The script is caltech_results/IA121_T9L62/run_arm.sh; the full key line is in each run.log.
- No rst was written. Each arm still wrote one end-of-run bin set (5.8 GB, kept in t9arms_1123/<arm>/bin, for comparison).

## Result, cycles 6450-6850
| | S (floor_pay_skip=false) | L (dfloor_flux_limit=false) | production IA12.1 (both true), same cycles |
|---|---|---|---|
| CLEAN WARNING | **0** | 38 | 80 |
| stage fallbacks | **0** | 18 | many |
| DT_RETRY | **0** | 22 (c6470 x7, 6530 x5, 6580 x3, 6625 x4, 6734 x1, 2 recovered) | 13 (6469 x3, 6498 x2, 6577 x8) |
| implicit_pos_floor created / energy unpaid | 0 / 0 | 1 / 1 | yes / 0 |
| dt min | **1.094** | 0.0278 | 0.0130 |
| dt at c6850 | 9.40 | 6.99 | ~7.6 (c6800) |
| t at c6850 | 4.77299e5 | 4.76105e5 | |
| CLEAN FALLBACK lines | 393 | 393 | |
| s/cycle (c6500-6850) | 3.09 | 4.67 | |
| Lbot/Lin / Ltop/Lin end | 0.268 / 2.532 | 0.396 / 2.546 | 0.373 / 2.544 (c6800) |
| max abs(skip_de) per [update] line | 2.3e42 (counted, not paid) | 4.6e42 | |

- Reading, by your table: **S clean and L like the production run, so the cause is paying the large floor energies from E_rad** (hydro/floor_pay_skip=true).
- Note that S's Lbot/Lin falls further (0.27) than L's or production's. This is the inner wall drift seen in every t9l62 arm from this state.

## Files: caltech_results/IA121_T9L62/arm_{S,L}/
- cw_dtretry_fallback.txt.gz: all CLEAN WARNING / DT_RETRY / CLEAN FALLBACK lines;
- floor_audit_pay.txt.gz: floor_audit [update] / [refresh] and the pay lines;
- clean_cycle_t2s.txt.gz: CLEAN cycle t2s lines;
- t9l62.user.hst.gz, t9l62.hydro.hst.gz, log.0.gz, run.log.

Node B (and node A, D0 done) are free again. Production t9l62 IA12.1 (.39), m35mms (.25) and the IA11 rst were not touched.

-- caltech
