# He giant scout128 (128^2, 30 d) - viper apu chain (user GO 10-07)

Setup and physics: /viper/ptmp2/jinma/he_giant_1006/SETUP.md.

- Binary: /viper/ptmp2/jinma/builds/bin/athena_he_gpu72_f7291255_hegiant, md5 ae8a3d11e0977dfe96f7fec7bd4bfbd6
  (branch hegiant-1006 f7291255, worktree /viper/ptmp2/jinma/wt_hegiant: he_gm_column + he_base_heat).
- Input: hegiant_scout128.athinput (md5 00e7dfbc). Differences from runs/inputs/hegiant_scout128.athinput: outputs
  (hst 25 s, hydro_w/m1/m1_vet bins 0.25 d, rst 0.5 d, log 1e4 s) and meshblock 445 x 64 x 64 (4 blocks).
  32 x 32 blocks FATAL in vet_gd (clamped lateral reads, 1.7e5 at cycle 0; smoke 12108597); 64 x 64: 0 clamped reads.
- run.cfg: BIN, BIN_MD5, IN, TLIM 2.592e6 (30 d), XKEYS "time/cfl_number=0.3" (named in IN), PK empty.
- Layout: 4 blocks, 1 per GPU, apu 2 nodes x 2 MI300A (8 GPUs cannot be used with 4 blocks).
- Seed: he_seed 1e-2, 16 modes, k 1-2, 3.2-55 Rsun (seed rule). Scaffold frozen, ramp 5 -> 10 d, base heat 0.5 Rsun.
- Chain script: dual_hegiant.sh (copy of /viper/ptmp2/jinma/bsg_viper_1006/dual_viper.sh made by mk_dual.py;
  only change: the first link starts fresh from IN when rst/ is empty and no run log exists). Guard, CANCEL,
  DONE, STOP, md5 check, linkcheck.py and rst_info.py from bsg_viper_1006 (read only).
- Smoke (rule 09-30): job 12108649, apudev 1 node, the same script/binary/input/run.cfg, all 4 blocks (2 per GPU),
  30 cycles: rc 0, FATAL 0, NON-CONV 0, Picard mean 4.2 max 7, dt 20.05 s constant, 1.26 s/cycle, band-clamped
  reads 0. Dir ../smoke128b.
- Jobs (jobs.txt): 12108729 (link 1, fresh) -> 12108730 .. 12108734, apu --nodes=2 --time=06:00:00, afterany.

How to stop / change:
- stop after the current link: `touch /viper/ptmp2/jinma/he_giant_1006/runs/scout128/STOP` (pending links exit 1);
  `CANCEL` = pending links exit 0. To kill the running link: scancel it on purpose (check its state first).
- change keys for pending links: edit run.cfg (XKEYS), then re-smoke on apudev before the next link starts.
- Quick look: `python3 /viper/ptmp2/jinma/he_giant_1006/runs/scout128/check_link.py`.

## Link 1 check (12108729, 2 nodes, 5.7 h wall incl. startup; checked once)
rc 0, FATAL 0, NON-CONV 0, Picard mean 4.3 max 19; 23400 cycles to t 4.467e5 s (5.17 d), 0.862 s/cycle on
4 MI300A, dt 20.0 -> min 16.6 -> 18.3 s. Newest rst hegiant.00011.rst. R_ph (shell mean, tau_R 2/3 with T_rad)
61.91 -> 62.81 (1 d) -> 64.44 (2.5 d) -> 65.87 Rsun (5.17 d), about +0.8 Rsun/d (~9 km/s); r(rho 1e-12) 62.97 ->
69.4 (3.5 d) -> 67.4 Rsun; max shell-mean v_r 40 km/s at 1.5 d (65-66 Rsun, the outer atmosphere), 20 km/s at 5.17 d.
Projection: 24.8 d left at dt ~18 s and 0.86 s/cycle ~ 33 h = ~6 links of ~5.4 h code time; 5 links are queued
(12108730-34), so one more link may be needed near the end.
