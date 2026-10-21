# NOTE 2026-10-21 caltech -> raven: accretor v12lzc_be_gm FRESH STARTED on 077dda23 (race with run12)

- Link 1 5223354 started 10-21 ~06:0x PDT on hpc-sm-02-04 (2 H200); chain body START 06:34:25 PDT after the in-job check: C2 vs C and R vs C ALL IDENT (all hst + bins).
- Static phase to t=0.05 done 07:05 PDT (4447 cycles, 1874 s, rc 0); stream on from static/rst/ryper.00005.rst.
- Stream: ~1.7 cycles/s (~0.6 s/cycle) at cycle ~8000; 0 NON-CONVERGED/FATAL; GPU peak 15.8 GB.
- Binary 077dda23 md5 63bca92a, deck 6b771411 (template 5de0635e). Links 5223355-59 afterok. Dir /resnick/groups/carnegie_poc/jingze/acc3d_v12_1021.
- Same-env question (TASK 1ba37054): this chain is 1 node, so the GDR setting doesn't matter for it; UCX_RNDV_THRESH=inf is ~4x slower on H200 (rdiff_ucx variant a), so not applied.
