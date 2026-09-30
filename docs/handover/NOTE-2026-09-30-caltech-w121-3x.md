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

**Update 09-29 17:35 PDT (user):** rebuilt at **793e03c3** (cubed-sphere seam-flux MPI request leak fix,
NOTE-2026-09-29-deltaai-seam-mpi-leak.md; waits only, bitwise unchanged): build 3643011 overwrites `athena.gpu`,
production 3640302 now afterok:3643011. The Caltech binary therefore differs from viper's 2fd94098 build by the
seam-wait fix only.

**Update 09-29 18:05 PDT: Caltech 3x STARTED FIRST.** Job 3640302 running since 09-29 17:46 PDT on hpc-sm-02-03
(binary 793e03c3). In-job smoke `SMOKE 3x rc=0 fatal=0 notconv=0`. dt 15.19 s at cycle 0 (= viper smoke), settles to
~12.0 s (11.6-12.4) from cycle ~4000 onward; ~16 ms/cycle, ~24.6 rot/h -> rot 300 ETA ~09-30 06:00 PDT (12.2 h, fits the
15:45 wall; no chain needed). No NaN/FATAL at cycle 56000 (rot ~6). Viper 12030356 had not started at this time ->
per the TASK the user cancels viper 12030356 (+12030357).

**Update 09-30 10:45 PDT: link 1 done, link 2 queued.** dt fell 12 -> ~7.5 s after rot ~18 (radial CFL at p ~1e-5 bar,
supersonic nightside downflows; see ck A/B TASK: not a solver limit cycle). Link 1 3640302 COMPLETED at wall (15:45,
rc 0) at t = 2.8017e7 s = rot 254.3 (cycle 3.517M, dt 7.47 s, no NaN/FATAL, last rst dhj.00509.rst, 2.9e7 zone-cyc/s).
Chain link 3648795 (afterok, same binary/keys, restarts from the newest rst) pending, est start 09-30 16:49 PDT;
~46 rot left (~3 h).
