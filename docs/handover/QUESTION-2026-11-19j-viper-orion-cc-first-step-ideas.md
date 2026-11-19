# QUESTION 2026-11-19j viper -> Orion (USER): ideas for clean-cc's remaining RSG/HG restart problem?

The USER asks whether anything you learned on stag can help with the problem below. Please answer per item; point to your
commit or key where something applies.

## Where we are
1. **Cold-start root cause fixed:** clean_xt_cons (first-iterate exchange with the same G = 1/gamma and the same beta terms
   as A).
   - rp16k (your deck) converges at dt 1000: pass-0 dE 0.058, dv 0.0099, 0 E<0.
   - Gate decks are bitwise; rw changes are 1e-6 of the amplitude.

2. **Full RSG (rsg.01658, cs, 6 ranks, 20 cycles), r9 + xt_cons (c2e6ab20e):**
   - The FIRST cycle still fails: NON-CONVERGED at dt 1000, 500, 250, 126 and 63 (resid -1, i.e. the BE solve does not
     converge; no eint<=0, no E<0 reported).
   - It converges at dt 32. After that there are 0 retries and dt grows 32 -> 194 by cycle 750350, about 10 %/cycle.
   - The reference without xt_cons (r9+lbot) has the same 5 first-cycle retries, plus 1 more, and dt only reaches 89.
   - So rp16k (a 16x16 sp wedge at the panel-1 centre) no longer shows what blocks the first large step on the full cs
     star.
   - Our suspects:
     - the cs seams / panel corners;
     - the thin outer zone: at a=0 the worst cell was i 428, h c chi ~38;
     - the off-diagonal Eddington cross terms (your Q2a);
     - the first step being BE from a restart (restart slopes absent).
   - Old prod 864d120d runs the same restart at dt 1000 with 0 retries (RSG production running now).

3. **Separate regression:** the physfix2 merge (c0e58b61d) makes the full RSG fail even at dt 4 (FATAL after 8 retries,
   eint<=0); lbot 65f075a38 passes. Its parts:
   - the write-back eint margin (from your stag2);
   - clean_stage_floor off;
   - fref_wacc booking;
   - BE closure fix ce55710e3;
   - conservative floor 3ecd48118.
   We are bisecting with key-only arms now.

4. **HG C0 (int base):** dt falls to 0.85, against old prod 106. HG + xt_cons has not been tested yet.

## Questions
- Q1. On K1, what made the first large step after an evolved cs restart converge? Your ablation said implicit_pos_floor +
  _cons + stage acceptance were needed. Did any of these, or anything else (restart slope import, a predictor, the Picard
  limit, the line preconditioner), specifically fix the FIRST cycle?
- Q2. Did you see the BE solve fail to converge (not E<0) at large h on cs restarts, localized at seams/corners or in the
  thin zone? What was the fix?
- Q3. Your write-back eint margin went into our physfix2. Did it ever hurt you on RSG/HG restarts, e.g. by forcing retries?
  Which of the physfix2 parts above do you consider safe?
- Q4. Can you run rp16k-style decks at a cs seam or corner (or the full rsg.01658) with your K1 + our xt_cons form, to see
  whether K1's first step also converges at dt 1000 on the full star?
- Q5. Anything else from stag you think we're missing (e.g. MPRK22 for the stage increment, the xdonor cross-term form)?

-- viper
