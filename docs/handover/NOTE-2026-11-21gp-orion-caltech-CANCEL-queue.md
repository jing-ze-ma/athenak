# NOTE 2026-11-21gp orion -> caltech (USER via coordinator): CANCEL all not-yet-started arms -- only H200CHK (NOTE gn + go) runs

The USER: "don't do unnecessary tests." Every BSG / T19 arm on B_v203 carries the suspected H200 closure trap, so:
1. **Cancel every arm that has not started yet**, and do not start them:
   - BG_bsgold and BH_bsgnowb;
   - BI_bsgtau1000 (if not started), BJ_bsgtau30 and the conditional BK_bsgtau1000thin0;
   - all H2* (H2A_* reruns, H2S_hg chain, H2C_hg, H2K/H2M), the V204S/T/U HG chain, the rest of v204c (CM*, V204*), PCS_*;
   - RSGDIE (RB_xpat, RA_on2, LIA4/LIA8), EW pairs, floor_pay seeds, ST/STD_t9l62, WGI_t9l62_corona and everything else in paused.txt.
2. **DGL_t19cutr:** stop it now if it is still running. Commit what it has so far (hst / raycut / log tail).
3. **Running BB / BD / BC / BE / BF** may finish (≤ ~30 min left). Commit their results as planned, but **start nothing new**.
4. **The ONLY new work is H200CHK:** NOTE gn (tiny 1-GPU BSG smoke + B_v203 build provenance) and its pre-authorised follow-up **A or B** from NOTE go. Commit to caltech_results/H200CHK/.
Please confirm the cancellation in one line: which arms were cancelled, which stopped, and which are still running.

-- orion (coordinator, for the USER)
