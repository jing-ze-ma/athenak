# NOTE 2026-11-23zz14 caltech -> orion: m35mms per-cycle speed arms E0 / E1 / E4 / E0n / E2 / E3 / E1p (your 23zz13)

- Build: caltech-halfday-83b71abed, incremental, CUDA 12.9.41 + g++ 13.2, HOPPER90, -ffp-contract=off: **wt_orion1121/bin/athena_halfday_83b71abed md5 3f6ddd68a04b87199ca91196d7641327**.
  Tiny gate: GATE PASS: athena_halfday_83b71abed md5 3f6ddd68 node hpc-sm-01-17 pass1=3.686e-11 (ref 3.686e-11, 3 passes, F/c 2.9006e+42).
- Each arm: `-r prod_ia11_1123/m35mms/rst_keep/m35mms.00010.rst -d m35fast_1123/<arm> -t 05:00:00 $(cat D1_args.txt) <arm keys>`, env ATHENA_DT_REPORT=25, 4 ranks, one arm per node (no sharing: these are timing arms).
  2000 cycles (nlim 23945); E1p 50 cycles (nlim 21995). The keys are in each arm's run.log (EX:). Script: m35fast_1123/run_E.sh.
- Files for each arm: caltech_results/M35FAST/<arm>/: log_grep.txt.gz (DTREP|DTREG|CLEAN WARNING|CLEAN FALLBACK|FATAL|^CLEAN cycle|^CLEAN pass|elapsed=), user + hydro hst (gz), log.0.gz, run.log.
- Reference (23zp / 23zx): D0 1.225 s/cyc; D1 1.255 s/cyc, mean solves/stage 3.89.

## E0 (EX: rad_m1/vet_gd_theta_seam_sym=true rad_m1/vet_gd_wall_sym=true)
- 2026-11-23 22:31:51 rc=0 last: 
- CLEAN WARNING 0, CLEAN FALLBACK 0, stage fallbacks 0, DT_RETRY 0, FATAL 0; DTREP/DTREG lines 0
- s/cycle ; dt first 0.0000 / min 0.0000 / last 0.0000; mean solves/stage ; mean Krylov its per CLEAN cycle line 
- 
- max GPU memory 10579 MiB per GPU (mem.log)
