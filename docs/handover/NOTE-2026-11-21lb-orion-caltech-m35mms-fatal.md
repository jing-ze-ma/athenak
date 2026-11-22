# NOTE 2026-11-21lb orion -> caltech: URGENT — m35mms FATAL (CW 33, dt 2.8, 200-pass stages); t9r CW 22, t9l62 CW 18

PROD 16:16: m35mms chunk ended with FATAL 1 after CW 33 and dt down to 2.8; t9r CW 3->22 (dt 13.6, 200-pass stages); t9l62 12->18.
NOTE ky (15:57) is not acked yet. Please do it now, in this order:
1. m35mms: restart from its newest GOOD rst (before the FATAL chunk) on B_p1 (branch caltech-prod-bf4f917d4, CUDA 12.9 + g++13,
   BSG smoke gate) with the kp keys + rad_m1/clean_gas_newton_switch=4 + clean_gas_newton_switch_min=6, env VGD_XFE=0, cfl 0.9.
   If it fails again within 100 cycles: same at cfl 0.6.
2. t9r, t9l62, f1k: same binary + keys (t9r cfl 0.6 + tfloor 1000 K; others 0.9).
3. Post m35mms's FATAL message + the 40 log lines before it -> caltech_results/PROD1121b/m35mms/fatal_1616.txt.
If building B_p1 blocks you (OTP / holder handover), say so in a NOTE right away.

-- orion (for the USER)
