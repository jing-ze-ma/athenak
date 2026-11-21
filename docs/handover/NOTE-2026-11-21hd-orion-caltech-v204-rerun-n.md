# NOTE 2026-11-21hd orion -> caltech (coordinator): rerun the v20.4 key A/B on B_v204n (CUDA 12.9), after BV + R_t19

The earlier v20.4 key evidence came from CUDA 12.2 binaries (B_v204c). That covers CM5_hg FATAL, CM4_hg with 8 CW and dt 6.8, CM5_rsg at +43 %, and V204_rsg / V204O. **Treat it as void.**
**After the BV arms and the R_t19 arms**, on test nodes, using **B_v204n** (caltech-runs/v204n 207fab564, CUDA 12.9, smoke-gated first), with **20 cycles each** and the same restarts as before (rsg.01295 / hegiant.00139; HG with time/dt_retry_host=m1):
- **CM4_hg**: hg_v204.keys (option A) + rad_m1/clean_conv_max=1e-4
- **CM5_hg**: hg_v204.keys (option A, conv_max 1e-5)
- **CM5_rsg**: rsg_v204.keys (conv_max 1e-5)
- **V204O_hg / V204O_rsg**: v20.3 keys + rad_m1/clean_maxit=200 rad_m1/clean_conv_max=1e-4. Expected: **bitwise** to a **B_v203n V203L rerun**. Please run V203L_hg / V203L_rsg on B_v203n for the same 20 cycles as the reference.
Report per arm, with the fixed parser: passes/cycle, s/cycle, CW, min dt, plus the hst compare for V204O vs V203L(n). These results decide the USER's v20.4 key choice.

-- orion (coordinator)
