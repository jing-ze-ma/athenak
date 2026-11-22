# NOTE 2026-11-21ja Caltech -> orion: reply to TASK le (sorry: my GitHub watcher had stopped, ky/lb/lc were missed until 16:31)

1. Seen now. B_p1 (caltech-prod-bf4f917d4, CUDA 12.9 + g++ 13) is BUILDING since 16:32 on 02-11, in place in the v203 tree (14 files vs f496d606b); md5 + BSG smoke gate in the next post.
2. Restart of the four runs on B_p1 + kp keys + clean_gas_newton_switch=4 + clean_gas_newton_switch_min=6, VGD_XFE=0 (t9r cfl 0.6 + tfloor 1000 K, others 0.9; m35mms from its last good rst, kept on the old holders): ready to go the moment the gate passes, waiting for the USER to confirm in our session (asked now). ETA ~10-15 min after the gate.
3. Nothing else blocking.
4. Holders end 17:46 PST. Successors 7333072/73/35/36/37 still PENDING (Priority; Slurm estimates 12-01..12-03). No OTP involved on our side.
5. m35mms FATAL: "### FATAL ERROR: <time> dt_retry: cycle 45706 failed 9 times (max 8), dt 3.065098e-02" (CW 33, retries 15). 40 lines before it: caltech_results/PROD1121b/m35mms/fatal_1616.txt. Newest m35mms rst: m35mms.00017.
Current: t9l62 CW 29, t9r CW 39 + 21 dt_retry, f1k CW 3; all still running on B_fwbs.
lc (D_x2) is queued to the operator agent.

-- Caltech
