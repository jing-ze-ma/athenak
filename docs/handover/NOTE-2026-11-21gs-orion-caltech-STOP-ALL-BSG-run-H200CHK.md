# NOTE 2026-11-21gs orion -> caltech (coordinator, USER: no unnecessary tests): STOP BB / BD / BC / BE / BF / BI NOW -- run H200CHK immediately

1. **Stop BB_bsgcutr, BD_bsgthinf0, BC_bsgtops, BE_bsgideal, BF_bsglq10 and BI_bsgtau1000 now** (step cancel). They all run on the suspected trapped B_v203 binary and are unnecessary per the USER. Commit their **partial hst** (plus raycut if handy); no further analysis is needed. Delete their own rst / bins.
2. **Run the H200CHK smoke + gq matrix immediately on the first freed node:** the tiny 1-GPU smoke of NOTE gn, run plain, with CUDA_LAUNCH_BLOCKING=1, and under compute-sanitizer memcheck, plus VGD_CS_LB=0 if it traps. Then the pre-authorised follow-up A or B (NOTE go). Commit to caltech_results/H200CHK/.
Nothing else runs.

-- orion (coordinator)
