# NOTE DeltaAI: BSG true repro R / RC smoke + cycle-0 IC check (TASK-2026-10-11-deltaai-bsgrec-smoke)

**Picked up** 2026-10-11 08:52 CDT.
- Built prodfinal-1011 **530c38e2** (= tip at build time): `athena_hes_gpu_530c38e2c428`, md5 **195a0a2fb73f509561b1e41df6b07014**.
- Files: our 10-09 files, md5 OK vs MD5SUMS; bsgR / bsgRC from the .in templates (diff = only implicit_gas_src + comments).
- ic0check.py sanity check on our own 10-09 GH200 hr cycle-0 dumps (bsg_hrdet_1009 job 3347471 detH_a, binary 98835d99):
  **BITWISE vs ic0_ref_1009.json in all 9 variables** (dens, velx, vely, velz, eint, m1_e, m1_f1, m1_f2, m1_f3).
- Job **3359915** (ghx4-interactive, 1 node x 4 GH200, 0:45): R4, RC4 (4 x 4, nlim 10), R2 (2 x 8, nlim 4), stop on rc != 0 / FATAL,
  nvidia-smi every 2 s. Extra command-line keys (both in the deck): `output1/dt=1e-3` (hst every cycle) and `output5/dt=1e30` (no rst).
  **PENDING** (priority 1987) at 09:00 CDT. STARTED line follows here when it runs.
