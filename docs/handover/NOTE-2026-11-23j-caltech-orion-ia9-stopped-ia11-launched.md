# NOTE 2026-11-23j caltech -> orion: IA9 stopped, IA11 built + gated + LAUNCHED (fresh bsg / m35mms / t9l62)

Reply to 23i. The USER confirmed "yes, stop IA9 and launch IA11" in session (relayed by COORD).

## 1. IA9 stopped (11-23 10:17 PST, clean)
`touch <run>/STOP`. Each run finished its cycle, wrote a final rst, printed "Terminating on stop request" and exited 0. Only the steps ended;
holder 7333072 was left alone. The IA9 poster (prodtab_ia9.sh) is stopped.
Final rows (the hst, chunks.txt and log.0.gz are in caltech_results/PROD_IA9_1123/<run>/, row in final_prod_row.txt):

| run | cycle | t | % tlim | dt | Ltop/Lin | Lbot/Lin | CW | fallbacks | DT_RETRY | MCAP lines | FATAL |
|---|---|---|---|---|---|---|---|---|---|---|---|
| bsg | 22500 | 3.431394e+06 | 69.68 % | 153.5 | 1.0237 | 1.0102 | 0 | 0 | 0 | 1280 | 0 |
| m35mms | 15850 | 1.274101e+05 | 14.75 % | 7.32 | 0.9880 | 1.0002 | 0 | 0 | 0 | 9683 | 0 |
| t9l62 | 6450 | 4.757564e+05 | 17.74 % | **4.96** | **2.5505** | **0.5059** | **2** | 0 | 0 | 12781 | 0 |

At the stop, IA9 t9l62 was going bad: dt had fallen from 96.8 to 4.96, L_top/L_in was 2.55 and L_bot/L_in 0.51.
Its two CLEAN WARNINGs are both floor_pay:

    CLEAN WARNING cycle 6381: floor_pay = rad_column left 2.096161e+40 of the fluid floors' energy unpaid (created)
    CLEAN WARNING cycle 6460: floor_pay = rad_column left 1.724604e+39 of the fluid floors' energy unpaid (created)

**IA9 rst: ALL KEPT, deletion waits for the USER's OK.** prod_ia9_1123/<run>/rst/:

| run | file | size | t | md5 |
|---|---|---|---|---|
| bsg | bsg3d.00000.rst | 21.8 GB | 0 | - |
| bsg | **bsg3d.00001.rst** (latest, final) | 22.3 GB | 3.431394e+06 | df39a7bbd946c2d6a5cb70dbf58875c2 |
| m35mms | m35mms.00000.rst | 28.5 GB | 0 | - |
| m35mms | **m35mms.00001.rst** (latest, final) | 29.1 GB | 1.274101e+05 | c27660bc9c04cc86f7d2523a179c2b38 |
| t9l62 | t9l62.00000.rst | 86.3 GB | 0 | - |
| t9l62 | **t9l62.00001.rst** (latest, final) | 88.1 GB | 4.757564e+05 | 6c5be78622405affc5c308abc916a0c8 |

**Unexpected, please note:** your IA9 and IA11 decks carry the CALSMOKE outputs. <output1> is hst every cycle; <output2> is rst with dt = 1e30
(final rst only). They have no bin outputs and no log output. So the IA9 runs wrote NO periodic rst and NO 3-D bin data in 4.8 h, and
00001 is their only end state. See section 3 for what I did for IA11.

## 2. IA11 binary + gate
- `caltech-ia11-bd9553537`, built in place (v203/build_n, incremental from 77551e823, 53 TUs, 228 s) with CUDA 12.9.41 + g++ 13.2.0,
  HOPPER90 and -ffp-contract=off. It ran as a CPU-only step in holder 7333072.
- `wt_orion1121/bin/athena_ia11_bd9553537`, **md5 bbfdee63ca6bf55897ea1886d0bc2391**.
- Tiny BSG gate (1 H200, hpc-sm-01-17, new check prod_ia11_1123/gate_ia11.sh; gate2.sh for IA9 left as is): **PASS**, identical to your A100
  result (caltech_results/H200CHK/gate_ia11.txt):

      pass 0 9.534e-02 | pass 1 3.686e-11 | pass 2 1.026e-14 -> converged after 3 passes; F/c -> 2.9006e+42

## 3. IA11 decks
Your caltech_ia11/IA11_*.athinput (md5 OK against MD5SUMS) with two changes:
1. The same file-path swap as IA9: he_ic_file, he_opac_table and he_planck_table (plus he_ic_linhao_file for bsg), pointed at our copies in
   fresh_prod_1120. The diff against our IA9 decks is exactly your 4 rad_m1 lines plus the header.
2. **Outputs: I replaced your CALSMOKE <output2> (rst dt 1e30) with the PROD1121b production outputs <output2..7>.** Those are bin hydro_w and m1
   (bsg 0.5 d, m35mms / t9l62 3 h), log 1e4, periodic rst (bsg 1 d, m35mms 6 h, t9l62 12 h), and m1_vet / m1_face bins (4x the bin dt).
   Your <output1> (hst every cycle) is kept. Physics, IC, <time> and <rad_m1> are unchanged.
   I did this per 23i sect. 5 ("same rst/bin cadences are fine"): a 3-day production with no periodic rst and no 3-D data seemed unintended.
   Say if you want it reverted.
- Decks as run: prod_ia11_1123/decks_final/, also in caltech_results/PROD_IA11_1123/decks/. md5 bsg dfa85cf9, m35mms d7e99b8f, t9l62 0b206728.
  Your deck with only the paths swapped is kept as orion_decks/IA11_<run>_caltechpaths_smokeout.athinput.
- nghost 3; meshblock nx1 = mesh nx1 in all three (no x1 split). NOTHING on the command line; cfl 0.3 comes from <time>.

## 4. IA11 launch (11-23 10:18 PST)
prod_ia11_1123/launch_ia11.sh STAR NODE 7333072 "2026-11-26 05:01" calls prodg.sh, i.e.

    OMPI_MCA_io=romio321 srun --overlap --jobid=7333072 -w <node> -N 1 -n 4 --mpi=pmix -J p_<run> \
      env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 memwrap.sh <dir> athena_ia11_bd9553537 \
      -i prod_ia11_1123/<run>/deck.athinput -d prod_ia11_1123/<run> -t <END - now - 10 min>

| run | node | step | layout |
|---|---|---|---|
| bsg | hpc-sm-01-17 | 7333072.12 | 4 ranks, 1 block/GPU |
| m35mms | hpc-sm-02-03 | 7333072.11 | 4 ranks, 1 block/GPU |
| t9l62 | hpc-sm-02-09 | 7333072.13 | 4 ranks, 2 blocks/GPU |

The holder ends 11-26 05:21. The expected wall is as in 23d: bsg ~6 h, m35mms 16-39 h, t9l62 19-30 h (more if dt falls).

## 5. First cycles (IA11, at 10:29) vs your 23i sect. 6
- 0 CW, 0 FATAL, 0 stage fallbacks, 0 DT_RETRY, 0 NON-CONVERGED in all three.
- 14 keys never read (the same set as IA9). The time2_tableau line prints no "time2_gtime_mom = false".
- user.hst has 30 columns (L_midE and L_intE added).

| run | cycle-0 dt | dt min (to c~750) | c200 dt | c200 Ltop/Lin | c200 Lbot/Lin | s/cycle | MCAP lines | 23i IA11 (A100) |
|---|---|---|---|---|---|---|---|---|
| bsg | 165.4553 | 164.64 | 168.72 | 1.0547 | **1.5341** | 0.775 | 2 | 165.46 / 164.64 / 168.72, 1.0547, 1.5341 |
| m35mms | 12.00764 | 12.0076 | 12.008 | 1.0011 | 0.9993 | 0.776 | 0 | 12.008, 1.0011, 0.9993 |
| t9l62 | 96.85234 | 96.850 | 96.851 | 1.0022 | **0.9387** | 2.32 (2 blocks/GPU) | 280 (to c~250) | 96.852 / 96.851, 1.0022, 0.9387 |

All shared rows match yours to the printed digits. clean_init converged in 4 passes for bsg (1.8e-13), and in 3 for m35mms (5.8e-12) and t9l62 (1.3e-9).
Outputs are working: bsg already has bsg3d.00001.rst (t = 8.64e4), and bin files are being written in all three.

## 6. PROD table
prod_ia11_1123/prodtab_ia11.sh appends a table every ~20 min to NOTE-2026-11-23k-caltech-orion-ia11-prod-table.md. Columns are as in 23g. The hst
files are copied to caltech_results/PROD_IA11_1123/<run>/.

-- caltech
