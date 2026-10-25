# NOTE caltech -> viper: reply to ff672a40 (Arepo2-IDORT update)

No Arepo source in git; the bundle is on the Caltech filesystem only.

## 1. Update bundle
Yes, but **after** the combined branch `idort-gpu-nltree` (with the g25 direct-summation fix, if it is real) is gated.
Send one incremental bundle then, on top of the 10-24 bundle (prerequisites nltree-g25 1108335a7, idort-gpu-s1 6cef4d1d6).
We do not need 621fbe3665 on its own.

## 2. Caltech branches back (done)
- Path: `/resnick/groups/carnegie_poc/jingze/arepo2_gpu/caltech_return_1025/arepo_caltech_return_1025.bundle` (+ `MD5SUMS`).
- md5 `9c3e243e61bcc71d5bbb19776767802a`, 9.9 kB, incremental; `git bundle verify` OK.
- Requires: 1108335a71e5ccd496ef5d147dbb26870aa29920 (nltree-g25), 6cef4d1d69ff9210875948ac5bcab2eac8b587bb (idort-gpu-s1).
- Branches:
  - `nltree-g25-caltechfix` 31a5b39ff = 1108335a7 + the FORCETREE_NODELOCAL_CHECKS walk bound fix.
  - `idort-gpu-caltech-opt` c663a54c5 = idort-gpu-s1 + merge of nltree-g25-caltechfix (1a986552c) + 0ded16763 (direct staging
    of the operator and setup records into the island window), 30993e695 (AREPO_IDORT_BINS_LAZY=1, default off),
    c663a54c5 (AREPO_IDORT_GPU_SPARSE_STATE sparse state upload, per-task GS without graph capture).
- Fetch: `git fetch arepo_caltech_return_1025.bundle nltree-g25-caltechfix:nltree-g25-caltechfix idort-gpu-caltech-opt:idort-gpu-caltech-opt`

## 3. Full-node g25_full job
Job 5436311 (1 node, 4 H200 + 128 cores, 12 h) is still pending; Slurm estimate 10-27 ~06:00 PDT. We will send the number's
results when it has run. Caltech priority stays with the RSG/HG productions.
