# NOTE 2026-11-21ga orion -> caltech (USER decision via coordinator): **STOP ALL production runs now -- BSG (01-17), T9L6.2 (02-11), m35mms (01-11)**

**The USER has decided to stop every remaining production run.** The USER may confirm in your session; please ask them there if your tooling needs it.
**Order of operations, per run:**
1. **First, before anything rotates or stops:** COPY (cp, never mv) the newest complete rst of each star, plus any pre-burst BSG rst still present (NOTE fz), into `/resnick/groups/carnegie_poc/jingze/prod_fresh_1121_saved/<star>/`. Verify each copy (size + md5 or cmp).
2. **Stop prodg.sh's KEEPRST pruning loop** for that star, so nothing more is deleted. Stopping the pruner is a stop, not a deletion.
3. **Stop the run cleanly.** Write a final rst if possible; otherwise cancel the step (as for T19).
4. **KEEP everything**: all existing rst files, logs, hst and bins. **No deletion of any kind.**

**Report per run** (one commit, final PROD_<star> refresh with full run.log.gz + hst):
- stop time and last cycle;
- every rst path kept (rst/ and prod_fresh_1121_saved/), with cycle and md5;
- final KE_int (hydro cols 8-10), L_bot/L_in and L_top/L_in (last user.hst row).

**Freed nodes** (01-17, 02-11, 01-11) join the test pool. Diagnosis / test queue first: PRODMON's T19 arms (DG, DE) and the BSG diagnosis arms (BSGDIAG, NOTE to follow). Then the rest per NOTE fg. One GPU job per node.

-- orion (coordinator, for the USER)
