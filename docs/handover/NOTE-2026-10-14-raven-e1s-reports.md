# NOTE Raven -> Caltech (10-14): E1S evolution reports + tlim rule applied

1. **E1S evolution reports (TASK-2026-10-14-caltech-e1s section 6) are needed from Caltech.** Raven cannot read the Caltech run dir,
   so nobody else can produce them. Please (whichever Caltech session can) run `pulse_ana.py run prod` + read user.hst and append one
   table per time to NOTE-2026-10-14-caltech-e1s.md as each is passed: t ~5e6, 7e6, 9e6, 1.25e7, 2.4e7 s, then every ~1e7 s.
   Columns as in section 6 (vrV, vrM, vr105, KEe, rsurf, L_top/L0, KE_int, w_mlt, Picard mean, DT_RETRY count, s/cycle, dt), next to the
   PU S3_268 / old-E1 reference rows. Push prod_pulse.txt + agcar3d.user.hst to docs/handover/e1s-1014/caltech_out/.
   Flag at once: vrV or vrM < -5 km/s, rsurf > 1.15, DT_RETRY bursts, Picard mean > 20.
2. **tlim rule (b3394647) applied on Raven/Viper:** He giant hegR (Raven 31135138-42) and hegV (Viper 12184016-18) run.cfg TLIM 2.592e6 -> 1e9
   (backups run.cfg.pre_tlim1e9; both chains still pending, no DONE files). Walltime there is already the max (24 h Raven gpu, 1 d Viper apu).
   Linhao BSG keeps 4.96e6 everywhere.
3. He giant page 9g8xcCyBCjdzwkBbALaiqt received, thanks. Mdot_top = 0: Raven is checking whether that column is wired for the He giant.
