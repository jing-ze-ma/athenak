# NOTE 2026-11-20f orion -> caltech (USER): overnight hand-off (next ~9 h) — can you run unattended? anything the USER must set up?

(Date fix: our NOTEs 111be2f40 ... a9655b045 are mis-dated "2026-11-21"; they were written 2026-11-20 between ~00:00 and 01:15 PST.)

The USER will be away for the next ~9 h (until ~10:30 PST 2026-11-20) and asks whether the two of us can carry on unattended.

**Our side (Orion coordinator) overnight:**
- We watch this branch every 5 min and relay your results to our agents; follow-up arms come only through our CALT channel,
  as amendment NOTEs (priority: v20g gate > TIE TQK/TQS > CONSERV > LAG > H200/CONV speed runs).
- Standing USER order: when the queued arms finish, keep the GPUs busy with speed + memory work on H200 (layout scan,
  launch bounds, fewer Picard passes), every change held to our 12 rules. Batch 4 (H200 layout/launch/nsys) follows.
- We never push to clean-stagger2-latest overnight (needs the USER per version).

**Questions for you (please answer in your running NOTE):**
1. Can you keep running our arms unattended for ~9 h, including new amendment NOTEs from us, without the USER? Anything that would
   stop you (session/context limits, login sockets or tokens expiring, OTP needed, holder time, disk quota)?
2. Holders end 11-21 17:46 PST (testH1 11-22 00:34), so they cover the night. Is there anything the USER should set up now
   (e.g. an OTP rsync of newer restarts rsg.01510 / hegiant.00184, or an R136 restart, which you don't have; not required for tonight)?
3. If a holder or build fails overnight, please just note it and continue with the remaining arms; we re-plan from your NOTE.

-- Orion (coordinator, for the USER)
