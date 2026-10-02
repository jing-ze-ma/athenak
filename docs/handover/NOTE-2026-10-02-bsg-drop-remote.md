# NOTE 2026-10-02 06:30 (viper): DROP the BSG runs on Caltech, DeltaAI, FREYA (and Raven) - user decision

The BSG arm-2 (and arm-1) productions STARTED on viper at 06:20 CEST (jobs 12059940 arm 2, 12059946 arm 1;
binary 30bf6c03; 1.1 s/cycle). The user decided to drop the BSG on all other sites.

Please, on each site:
- **Cancel** every queued or running BSG job (production and gate jobs) and do not create the READY marker.
  Caltech: production 3746142 and gate 3750607 (if still queued). DeltaAI: the queued 24 h ghx4 job and any
  interactive links. FREYA: do not start TASK-2026-10-02-freya-bsg (cancelled). Raven: handled by viper.
- **Keep** the builds, the gate outputs and the scripts (useful for a later cross-check); delete nothing else.
- Push a one-line NOTE confirming what was cancelled (job ids).
TASK-2026-10-01-caltech-bsg, TASK-2026-10-01-deltaai-bsg, NOTE-2026-10-02-bsg-queue-now and TASK-2026-10-02-freya-bsg
are closed by this NOTE. The comparison plan in NOTE-2026-10-02-sync-bsg is void.
