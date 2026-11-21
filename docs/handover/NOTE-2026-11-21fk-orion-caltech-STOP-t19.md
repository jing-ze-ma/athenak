# NOTE 2026-11-21fk orion -> caltech (USER decision via coordinator): **STOP T19L6.4 (02-09) and the T19 corona test (02-04) now**

The USER has decided to stop both runs: **t19l64** (hpc-sm-02-09, holder 6133591) and **t19l64_corona** (hpc-sm-02-04, holder 6133592). The reason is a numerical radial instability, per PRODMON's evidence.
- At 17:58: T19 KE_int 2.8e43 and doubling every hour, L_bot/L_in swinging between -3.2 and 2.2, dt falling from 38 to 24, and 136-168 eint<0 OLDVEC fixes per cycle. T19 corona KE_int 1.2e43 and L_top/L_in 15.
- **Your side will probably need the USER's confirmation in your session again.** The USER has been told and expects your prompt.

How to stop:
- Stop each run cleanly. Prefer a final restart dump: signal prodg.sh or the step so that AthenaK writes its final rst; otherwise cancel the step.
- **Keep everything**: the latest rst file(s), the full run.log (commit it gzipped as `PROD_<star>/run.log.gz`), and the full hydro/user hst (final PROD_<star> refresh).
- Commit one line per run: stop time, last cycle, the rst path(s) kept, and the final KE_int and L_bot/L_in.

**Keep T9L6.2 (02-11), m35mms (01-11) and BSG (01-17) running.** Do not touch them.

The 2 freed nodes (02-09, 02-04) go to:
1. T19 diagnosis / test arms, owned by PRODMON. They will come in a separate NOTE.
2. The stretched-grid T19 arms ST_t19l64 / STD_t19l64 (NOTE er, owner CORONA9), so move those up.
3. Then the queued tests in the NOTE fg order (v20.4c, ST/STD, WGI_t9l62_corona, CONVERR, H2*/PCS, EW, seeds).
The NOTE fj profiles stay wanted, from the kept T19 hydro_w bins.

-- orion (coordinator, for the USER)
