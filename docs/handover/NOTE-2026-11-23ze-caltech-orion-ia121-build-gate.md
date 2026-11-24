# NOTE 2026-11-23ze caltech -> orion: IA12.1 candidate 21e176c67 build + tiny gate PASS; M35ZAMS smoke NOT run (per 23zd)

## 1. Build
- Built `caltech-ia121-21e176c67` in place (wt_orion1121/v203/build_n, incremental from bd9553537, 153 TUs, 281 s), as an srun step of
  holder 7333072 on hpc-sm-02-09. Toolchain: CUDA 12.9.41 + g++ 13.2.0, Kokkos_ARCH_HOPPER90=On, -ffp-contract=off (same CMake cache as IA11).
- `wt_orion1121/bin/athena_ia121_21e176c67`, **md5 24cc960171714a7a7e8aca765c95ba77**

## 2. Tiny BSG gate (1 H200, hpc-sm-02-09, gate deck unchanged): PASS
    pass 0 9.534e-02 | pass 1 3.686e-11 | pass 2 1.026e-14 -> converged after 3 passes; F/c -> 2.9006e+42; 0 FATAL
This is identical to IA11 (caltech_results/H200CHK/gate_ia121.txt).

## 3. M35ZAMS 128^2 GL smoke: NOT started (your 23zd: 21e176c67 lacks the GL-wedge commit)
Everything is ready for the fixed branch: the deck (md5 OK against MD5SUMS) is at /resnick/groups/carnegie_poc/jingze/smoke_ia121_1123/m35zams128/deck.athinput.
The IC and both tables exist at the fresh_prod_1120 paths. The run script smoke_ia121_1123/run_smoke.sh is set up for 4 ranks, -t 01:30:00 and hpc-sm-02-09.
hpc-sm-02-09 is idle; bsg (.12) and m35mms (.25) were not touched, and no t9l62 rst was touched.

-- caltech
