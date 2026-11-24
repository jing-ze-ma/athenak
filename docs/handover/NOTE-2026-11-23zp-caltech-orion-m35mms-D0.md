# NOTE 2026-11-23zp caltech -> orion: m35mms arm D0 (your 23zm) done

- Ran 18:24-19:05 PST on N2 7333073 node A (hpc-sm-02-11, 4 H200). BIN athena_ia121_fa59bffd6 (md5 9b4e69b9); D0_args.txt md5 c8cf1535 OK; env ATHENA_DT_REPORT=25.
- `-r prod_ia11_1123/m35mms/rst_keep/m35mms.00010.rst` (cycle 21945, t 1.728e5; read only) `-d m35fast_1123/D0 -t 05:00:00 $(cat D0_args.txt)`.
- **The 3 mesh/ keys were accepted; no FATAL on layout_invariant_mode=tile.** layout_tile was reported "added from the command line".
- Ran 2000 cycles to nlim 23945, rc 0. **0 CLEAN WARNING, 0 CLEAN FALLBACK, 0 stage fallbacks, 0 DT_RETRY, 0 FATAL**; 160 DTREP/DTREG lines.
- **s/cycle 1.225** (cycles 21950-23900). For comparison, the production IA11 run took 1.134 s/cycle over the same cycle range (prod_ia11_1123/m35mms log).
- dt 8.61 (c21950) -> 7.19 (c23900); at c23900 production IA11 has 7.193585 and D0 7.193016. Mean 4.19 solves per stage.
- Ltop/Lin 0.98-1.004 and Lbot/Lin 0.996-1.007 (every 400 rows in the hst).
- GPU memory per GPU: see mem.log in the run dir (sampled at 2 s).
- Files in caltech_results/M35FAST/D0/:
  - log_grep.txt.gz: `grep -E "DTREP|DTREG|CLEAN WARNING|CLEAN FALLBACK|FATAL|^CLEAN cycle|elapsed="`;
  - m35mms.user.hst.gz, m35mms.hydro.hst.gz, log.0.gz, run.log, run_D0.sh.
  The bins (9.8 GB) are kept in m35fast_1123/D0/bin; no rst was written.
- Node A is free again. Node B is running the t9l62 arms S, then L (23zn).

-- caltech
