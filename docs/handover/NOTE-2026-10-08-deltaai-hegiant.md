# NOTE 2026-10-08 (DeltaAI): He giant N897 per TASK-2026-10-07-deltaai-hegiant -- smoke PASS, link 1 queued

- **Bundle.** `hegiant_deltaai_bundle.tar` md5 6acf7ec374bb90fb98deb18f46dbe6e1 OK. It is unpacked to
  /work/nvme/bivj/jma20/hegiant_deltaai_bundle, and SETUP.sh reports "md5: all 21 files OK".
- **Binary.** `athena_hes_gpu_f3a66907dc9c` (f3a66907, PROBLEM=he_star_m1, build_inc_deltaai.sh target hes_gpu),
  md5 **726ee6cd2435aa9c0c66539680055887**.
  - Stack: `module reset` default = PrgEnv-gnu 8.6.0, gcc-native/14, cudatoolkit/25.5_12.9, cray-mpich/9.0.1,
    craype-accel-nvidia90.
  - Flags: Kokkos HOPPER90 + ARMV9_GRACE, nvcc_wrapper, MPI on. FMA contraction is on (the default).
- **Smoke.** Job **3338065**, ghx4-interactive, 1 node (gh066) x 4 GH200, 4 ranks, 1 block per GPU, 60 cycles.
  The run dir is /work/nvme/bivj/jma20/hegiant_1008/smk60; its run.cfg, job output, check_link output and thinned hst are in
  `hegiant-deltaai-1008/`.

| quantity | DeltaAI 3338065 | viper 12120235 |
|---|---|---|
| rc / FATAL / NaN | 0 / 0 / 0 | 0 / 0 / 0 |
| dt at 47861 | 9.223353e+00 | 9.223353e+00 |
| dt at 47862 / 47871 / 47881 / 47921 | 1.844671e+01 / 1.856086e+01 / 1.856955e+01 / 1.858862e+01 | same digits |
| Picard (solves, mean, max, NON-CONV) | 60, 18.5, 42, 0 | 60, 18.5, 42, 0 |
| NEWTON-FALLBACK lines | 1920 | 960 (4 ranks here vs 2 on viper; informational) |
| last hst t | 9.0152570477933926e+05 | 9.0152570477933483e+05 |
| last hst dt | 1.8588619517720179e+01 | 1.8588619517698017e+01 (1.2e-12 rel) |
| last hst mass | 1.5435375539788512e+32 | 1.5435375539788510e+32 (1.3e-16 rel) |
| last hst tot-E | 1.5044198258655806e+47 | 1.5044198258655804e+47 (1.3e-16 rel) |
| he_ic_balance | 0.00515103 at r = 2.15997e+11 | same |
| R_ph(tau_R 2/3) / r(rho 1e-12) / max abs(<v_r>) | 73.316 / 76.048 Rsun / 37.579 km/s at 72.79 | same |
| s/cycle | **0.658** (median, cycles 10-60; check_link 0.701) on 4 GH200, 1 block per GPU | 2.779 on 2 MI300A, 2 blocks per GPU |

**PASS** on every criterion of TASK section 5.

- **Production link 1.** Job **3337674**: ghx4, 1 node x 4 GH200, 48 h, `-J hegiant897`, submitted 2026-10-07 23:37 CDT
  (before the smoke, with the user's OK, to hold queue position). Status: PENDING.
  - Run dir /work/nvme/bivj/jma20/hegiant_1008/run897. Its run.cfg equals the smoke's, with `time/restart_refill_ghosts=true`
    for link 1.
  - **Cost estimate:** about 91k cycles x 0.66 s = about 17 h of compute, so one 48 h link should reach TLIM 2.592e6 s,
    about 70 SU.
  - **After link 1:** drop refill_ghosts and re-smoke 10 cycles before any link 2 (TASK section 6).
- **Viper.** The viper chain hegiant897 is untouched; the user decides which copy continues.
