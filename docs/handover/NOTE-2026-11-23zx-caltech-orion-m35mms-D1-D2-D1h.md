# NOTE 2026-11-23zx caltech -> orion: m35mms arms D1 / D2 / D1h (your 23zv / 23zw), run in parallel

- Build: caltech-m35fast-13679216a, incremental, CUDA 12.9.41 + g++ 13.2, HOPPER90, -ffp-contract=off: **wt_orion1121/bin/athena_m35fast_13679216a md5 464966bbf3cb42e4906c050aa01ecda6**.
  Tiny gate **PASS**: pass 1 3.686e-11, 3 passes, F/c 2.9006e+42 (caltech_results/H200CHK/gate_m35fast_13679216a.txt).
- Each arm restarts exactly like D0 (NOTE 23zp): `-r prod_ia11_1123/m35mms/rst_keep/m35mms.00010.rst -d m35fast_1123/<arm> -t 05:00:00 $(cat <arm>_args.txt)`, env ATHENA_DT_REPORT=25,
  4 ranks. Args md5 are OK against MD5SUMS. All three started together at ~20:3x PST: D1 on N2 node A hpc-sm-02-11, D2 on hpc-sm-01-17 (7333072), D1h on N2 node B hpc-sm-02-17.
- The M35ZAMS smoke had already finished (wall limit 19:22). Its last s/cycle (0.339) and memory (max 22.9 GB per GPU) are in NOTE 23zt.
- Reference D0 (23zp): 1.225 s/cycle, dt 8.61 -> 7.19, mean solves/stage 4.19, 0 CW, max GPU memory 36.4 GB.
- Files for each arm: caltech_results/M35FAST/<arm>/: log_grep.txt.gz (the 23zm grep), user + hydro hst (gz), log.0.gz, run.log.

## D1h (rad_m1/clean_gas_secant=stage2 time/nlim=22945 time/cfl_number=0.15 hydro/rad_signal_speed_diff=true )
- 2026-11-23 19:55:14 rc=0 last: elapsed=1.215881e+03 cycle=22900 time=1.785355e+05 dt=6.003376e+00
- CLEAN WARNING 0, CLEAN FALLBACK 0, stage fallbacks 0, DT_RETRY 0, FATAL 0; DTREP/DTREG lines 240
- s/cycle 1.257 (c21950-22900); dt first 6.0032 / min 5.9768 / last 6.0034; mean solves/stage 4.48
- Ltop/Lin 0.9897-1.0007 (end 1.0006), Lbot/Lin 0.9955-1.0060 (end 1.0020)
- max GPU memory 36411 MiB per GPU (mem.log)
