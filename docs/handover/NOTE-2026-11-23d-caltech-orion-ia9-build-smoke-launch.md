# NOTE 2026-11-23d caltech -> orion: IA9 build + smoke + fresh bsg/m35mms/t9l62 LAUNCHED (per 23f GO)

Launch approved by the USER (relayed by COORD, orion NOTE 23f). Started 11-23 05:27 PST.

## Binary
- `caltech-ia9-77551e823` built in place (v203/build_n, incremental from 462210d51, 53 TUs, 224 s), CUDA 12.9.41 + g++ 13.2.0,
  Kokkos_ARCH_HOPPER90=On, -ffp-contract=off.
- `wt_orion1121/bin/athena_ia9_77551e823`, **md5 bb0f2fbe7030a0cb19f0bd054cfa6de7**

## Smoke (tiny BSG gate, 1 H200, same deck as before; caltech_results/H200CHK/gate_ia9.txt)
    pass 0 9.534e-02 | pass 1 1.885e-06 | pass 2 9.371e-12 | pass 3 2.887e-15 -> converged after 4 passes
    E 4.1527620104e+45 -> 4.1602693712e+45, F/c -> 2.9006e+42 1.1804e+39 -2.0678e+33
- Not the CUDA-12.2 trap (pass 1 3.5e-8, F small): E matches the healthy value (4.1604e45 within 4e-5) and F streams out.
- Different from head 462210d51 (pass 1 4.934e-06, NOT converged after 50, F/c 1.355e42). IA9 converges to round-off in 4 passes, and F/c
  is 2.1x larger. Please compare on your side (gate2.sh prints "FAIL" only because pass 1 != the old 4.020e-06 reference).

## Decks
Your caltech_ia9/IA9_*.athinput (md5 OK against MD5SUMS). The only change: the 10 viper file paths are replaced by our copies of the same
files, as in the PROD1121b decks (he_ic_file, he_opac_table, he_planck_table; plus he_ic_linhao_file for bsg) under
/resnick/groups/carnegie_poc/jingze/fresh_prod_1120/. The decks as run are in caltech_results/PROD_IA9_1123/decks/ (md5 bsg 35ffa908,
m35mms 52da1397, t9l62 ea1da0e6). cfl 0.3 comes from <time>. Nothing is on the command line.
mesh/nghost = 3 in all three, and every deck has meshblock nx1 = mesh nx1 (no x1 line is split).

## Command lines (prodg.sh via prod_ia9_1123/launch_ia9.sh STAR NODE HOLDER END)
Each run: 1 node = 4 H200, 4 ranks (as our previous productions; bsg/m35mms 4 blocks = 1/GPU, t9l62 8 blocks = 2/GPU). Holder h1120_N3
7333072 (ends 11-26 05:21; wall -t 71:23:52, END 11-26 05:01).
    OMPI_MCA_io=romio321 srun --overlap --jobid=7333072 -w <node> -N 1 -n 4 --mpi=pmix -J p_<run> \
      env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 memwrap.sh <dir> athena_ia9_77551e823 \
      -i prod_ia9_1123/<run>/deck.athinput -d prod_ia9_1123/<run> -t 71:23:52
- bsg    : hpc-sm-01-17 (step 7333072.5)
- m35mms : hpc-sm-02-03 (step 7333072.3)
- t9l62  : hpc-sm-02-09 (step 7333072.4)
Run dirs: /resnick/groups/carnegie_poc/jingze/prod_ia9_1123/<run>/. The t=0 rst files are 21 / 27 / 81 GB; 49 TB free.

## First cycles (log.0 + user.hst at 05:38)
| run | cycles | CW | stage fallbacks | DT_RETRY | dt | s/cycle | L_top/L_in | L_bot/L_in | orion smoke (s/cyc, dt, Lt/Lin) |
|---|---|---|---|---|---|---|---|---|---|
| bsg    | 750 | 0 | 0 | 0 | 165.46 -> 169.26 | 0.754 | 1.0227 (1.04-1.06 at c200-500) | 1.047 | 1.41 (2 ranks), 165.46, 1.0556 |
| m35mms | 700 | 0 | 0 | 0 | 12.0076 | 0.798 | 1.0003 | 1.0005 | 0.83, 12.008, 1.0021 |
| t9l62  | 220 | 0 | 0 | 0 | 96.852 -> 96.844 | 2.45 (4 GPUs, 2 blocks each) | 1.0024 | 0.8625 | 1.19 (8 ranks), 96.85, 1.0023 |
- Cycle-0 dt matches yours exactly in all three.
- clean_init: bsg converged after 5 passes (1.1e-13), m35mms after 4 (3.7e-12), t9l62 after 3 (1.4e-9).
- No FATAL. Startup warnings are the same as before: the arad / temp_unit "inconsistent with <units>" lines, and 14 keys never read.
  `rad_m1/clean_stag_thinf_held` is never read because thinf = false; the others are the usual cs/time/problem keys.
- MCAP: bsg 1 line (cycle 4, 228 terms); m35mms 0. t9l62 has the start transient you expected: 17k-60k terms per cycle at cycles 0-15, 1424 at
  cycle 17. After that it does NOT go away completely: from cycle 69 on, a few cycles per 15 print 8-128 terms (cycle 219: 48 + 48).
  Small but persistent. Please say whether that is expected.
- t9l62 L_bot/L_in is 0.86 (0.78-0.87 since t ~ 3e3 s); please check against yours.

## Expected wall (at the current dt; earlier productions' dt fell later, so these are lower bounds)
bsg ~29k cycles -> ~6 h; m35mms ~72k cycles -> ~16 h (J3 production pace: ~39 h); t9l62 ~28k cycles -> ~19 h (J3 pace: ~30 h).
All three fit in the 72 h holder.

## PROD table
A new poster (prod_ia9_1123/prodtab_ia9.sh) appends a table every ~20 min to NOTE-2026-11-23g-caltech-orion-ia9-prod-table.md. Columns:
CW, fallbacks, DT_RETRY, MCAP, FATAL, dt, s/cycle, L_top/L_in, L_bot/L_in and the last 3 hst rows. The hst files are copied to
caltech_results/PROD_IA9_1123/<run>/.

-- caltech
