# TASK raven -> Caltech: REQUEUE the REAL BSG true reproduction (user decision)

User: "we should requeue bsg on viper and caltech. delta ai is too slow".
The race-stop NOTE (NOTE-2026-10-12-caltech-racestop.md) is WITHDRAWN for Caltech.

Do exactly what TASK-2026-10-12-caltech-bsgprod.md says, unchanged:
- build bsgcand-1012 @ d3bf42d7 (report md5);
- deck docs/handover/bsgprod-1012/bsgP_seed.athinput.in (the 53995cf9 fix) with @FILES@ resolved only;
- smoke + ic0check, then fresh production to tlim 4.96e6 on 1 node x 2 H200 in chained links.

Copies now running or queued: DeltaAI 3367085 (running, keep it), Raven 31107715-17, Viper 12165594-96.
New rule: do NOT stop your copy when another site starts; all copies run until the coordinator says otherwise.
Report in NOTE-2026-10-12-caltech-bsg-requeue.md: QUEUED (job ids, est start), STARTED (md5, smoke, ic0check), then per link
(t, cycle, dt, s/cycle, Picard, counts, newest rst).
