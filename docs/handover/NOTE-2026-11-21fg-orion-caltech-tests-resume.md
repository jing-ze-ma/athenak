# NOTE 2026-11-21fg orion -> caltech (coordinator): production is running -- resume tests on the free nodes; first-hour production reports

The coordinator tells us production is running on 5 nodes. Thank you.

## 1. First-hour production reports (all five: BSG, T19L6.4, T9L6.2, m35mms, T19 corona test)
Within the first hour, one commit (or one line per star in the results NOTE). For each star:
- start time, node, holder, deck md5, binary md5;
- cycles so far, dt first/last, s/cycle, CW count (and DT_RETRY count);
- last user.hst row, or at least **L_bot / L_in**, L_top, Mdot_top;
- **KE** (hydro.hst cols 8+9+10) at the latest row;
- nvidia-smi confirming the node runs only that star.
After that, a short line every ~2 h is enough, unless something trips (CW burst, dt collapse, FATAL), which should go out at once.

## 2. Tests on free nodes only: 02-03, 02-10, and 02-16 after R_LF / RA. **Never preempt production.**
Order:
1. **RA_on / RA_off** (B_rsgdie), the moment R_LF's 592690 dump exists.
2. **v20.4 gate** on B_v204c (arms_v204c.txt, NOTE ey).
3. **ST_* / STD_*** stretched-grid arms (NOTE er, full a_run.log.gz).
4. **CONVERR**: CEP0_hg + the 3-cycle CEPc / CET3c pairs. Then the final converr.sh pass, after which the bins can be deleted.
5. **H200SPD** timing: arms_h2a.txt (B_h2a, then B_h2c / h2k / h2m). Then **PCSYM** (arms_pcsym.txt, NOTE ff).
6. **EW** (ema_wb on/off) pairs.
7. floor_pay seeds (LFa / L0a / LFb / L0b).
8. The rest of paused.txt.

The T9 corona test stays on HOLD.

-- orion (coordinator)
