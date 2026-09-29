# NOTE 2026-09-30 (Caltech): WASP-121b 3x duplicate production (TASK-2026-09-30-caltech-w121-3x)

**Start estimate (sbatch --test-only, 09-29 16:40 PDT):** 1 node x 2 H200 at **2026-09-30 15:23 PDT
(= 10-01 00:23 CEST)**. Viper's copy 12030356 was estimated at 09-30 07:35 CEST
(= 09-29 22:35 PDT), so **viper is expected to start ~17 h earlier**. Caltech jobs are queued anyway; if viper
starts first, cancel the Caltech ones (`scancel 3640301 3640302`, nothing spent while pending).

- Dir `/resnick/groups/carnegie_poc/jingze/w121prod_0930/` (package w121_3x_pkg md5 9dc9472ee4545b1b0150767577c068ec,
  SETUP.sh run with data from athenak_data/exo_fms_ck; input = the package w121prod_3x.athinput unchanged).
- Build 3640301 (expansion): build_inc.sh w3x gpu **2fd94098** -> `athena.gpu`, BUILD_OK, md5 0b507eeeb19ef2a6cd72ca832ffe685b.
- Production 3640302 (afterok build): 1 node x 2 H200, `prod.sub` ARM=3x (copy of the 1x/10x w121prod_0928 script),
  fresh start with a 50-cycle in-job smoke (exits 1 on FATAL / NOT-CONVERGED / rc != 0), then tlim 3.304605e7 (rot 300),
  -t 15:45. 10x took 13.7 h on 2 H200 here; if 3x does not reach rot 300 in one link, chain with afterok.
- Input diff vs the Caltech 10x input: only the expected 3x keys (grid nx1 76 / x1 / stretch coefficients, grav, ap, met,
  eos_xh/yhe, rad_met, ck tables, IC) plus ck_impl_conserve = 1, ck_impl_tol 1e-8, maxit 16, dtmax 0.5 (viper's input as shipped).
