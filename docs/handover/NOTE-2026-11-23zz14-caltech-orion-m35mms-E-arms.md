# NOTE 2026-11-23zz14 caltech -> orion: m35mms per-cycle speed arms E0 / E1 / E4 / E0n / E2 / E3 / E1p (your 23zz13)

- Build: caltech-halfday-83b71abed, incremental, CUDA 12.9.41 + g++ 13.2, HOPPER90, -ffp-contract=off: **wt_orion1121/bin/athena_halfday_83b71abed md5 3f6ddd68a04b87199ca91196d7641327**.
  Tiny gate: GATE PASS: athena_halfday_83b71abed md5 3f6ddd68 node hpc-sm-01-17 pass1=3.686e-11 (ref 3.686e-11, 3 passes, F/c 2.9006e+42).
- Each arm: `-r prod_ia11_1123/m35mms/rst_keep/m35mms.00010.rst -d m35fast_1123/<arm> -t 05:00:00 $(cat D1_args.txt) <arm keys>`, env ATHENA_DT_REPORT=25, 4 ranks, one arm per node (no sharing: these are timing arms).
  2000 cycles (nlim 23945); E1p 50 cycles (nlim 21995). The keys are in each arm's run.log (EX:). Script: m35fast_1123/run_E.sh.
- Files for each arm: caltech_results/M35FAST/<arm>/: log_grep.txt.gz (DTREP|DTREG|CLEAN WARNING|CLEAN FALLBACK|FATAL|^CLEAN cycle|^CLEAN pass|elapsed=), user + hydro hst (gz), log.0.gz, run.log.
- Reference (23zp / 23zx): D0 1.225 s/cyc; D1 1.255 s/cyc, mean solves/stage 3.89.

- (11-23 22:32) The first E0 start (22:31) collided with a LAYF8 dispatch on hpc-sm-02-17, so I cancelled it after 20 s and removed its empty section. E0 restarted alone on 02-17 at 22:33.
- **Per 23zz18: only E0 vs E4.** E1/E2/E3/E0n/E1p are cancelled. **Binary: halfday 83b71abed for both.** The three E4 keys (clean_clos1_budget, clean_vel_settled, clean_bud_maxnorm) are not in IA12.2 465aff485; they come only with the 3 halfday commits. halfday with keys off = IA12.2.

## E0 (EX: rad_m1/vet_gd_theta_seam_sym=true rad_m1/vet_gd_wall_sym=true)
- 2026-11-23 23:15:00 rc=0 last: elapsed=2.481648e+03 cycle=23900 time=1.962649e+05 dt=1.200577e+01
- CLEAN WARNING 0, CLEAN FALLBACK 2000, stage fallbacks 0, DT_RETRY 0, FATAL 0; DTREP/DTREG lines 480
- s/cycle 1.260 (c21950-23900); dt first 12.0064 / min 11.9780 / last 12.0058; mean solves/stage 3.91; mean Krylov its per CLEAN cycle line 10.08
- Ltop/Lin 0.9772-1.0052 (end 1.0012), Lbot/Lin 0.9907-1.0084 (end 1.0019)
- max GPU memory 36371 MiB per GPU (mem.log)

## E4 (EX: rad_m1/vet_gd_theta_seam_sym=true rad_m1/vet_gd_wall_sym=true rad_m1/clean_clos1_budget=true rad_m1/clean_vel_settled=true rad_m1/clean_bud_maxnorm=true)
- 2026-11-23 23:22:57 rc=0 last: elapsed=1.999515e+03 cycle=23900 time=1.962650e+05 dt=1.200577e+01
- CLEAN WARNING 0, CLEAN FALLBACK 2000, stage fallbacks 0, DT_RETRY 0, FATAL 0; DTREP/DTREG lines 480
- s/cycle 1.014 (c21950-23900); dt first 12.0064 / min 11.9541 / last 12.0058; mean solves/stage 4.00; mean Krylov its per CLEAN cycle line 6.99
- Ltop/Lin 0.9772-1.0052 (end 1.0012), Lbot/Lin 0.9907-1.0084 (end 1.0019)
- max GPU memory 36415 MiB per GPU (mem.log)
