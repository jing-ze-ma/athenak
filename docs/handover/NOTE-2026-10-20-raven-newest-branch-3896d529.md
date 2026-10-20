# NOTE 2026-10-20 raven -> caltech + deltaai: build and use origin/prodmerge-1019b @ 3896d529 NOW (USER: use the updated branch everywhere)

3896d529 = fin51-r-fix 0fbfa31f + torque registers (lz.hst) + Emsk booking fix + hr_q=col opt-in + time2_reject=retry BE fallback
(82275582) + opt-in <problem>/coriolis_lz (accretor). A100 gates clean (strict gbit 4/4 GMD5 0; accretor key on/off det/restart/retry
BITWISE); MI300A gbit pending; fin51-r-fix will be fast-forwarded to it after that. Raven + Viper chains already switched.
Why it matters: with 0fbfa31f/563ab474, time2_reject=retry can FATAL ("cycle N failed 9 times") on a dt-independent rejection;
3896d529 falls back to the BE redo after 2 cuts (<time>/dt_retry_t2_be, default 2; no deck change).

## Caltech
- RSG E copy (link 1 5126654 pending): rebuild he_star_m1 H200 nofma at 3896d529 and point the queued link at it BEFORE it starts
  (deck unchanged). Report md5.
- Accretor trc run (13206541, running): at its NEXT link switch to the 3896d529 accretor build (your nofma recipe) AND add one line to
  the deck under <problem>:  coriolis_lz = true
  (Raven's deck: v11amwtr + that line; inertial L_z about the accretor axis then closes to round-off, the Coriolis/transport drift goes away;
  energy: Emsk no longer books the sink KE, so mk_prod.py:90 is correct for hst written from then on). Mark the switch time on the page.
- The restart triplet 5140304 (563ab474): keep it, it still answers the H200 restart question.

## DeltaAI
- RSG E chain: rebuild at 3896d529 (GH200 nofma) and use it from the NEXT 2 h link on (restart from the newest rst; deck unchanged).
- Race rule for your copy (coordinator proposal, pending USER): a DeltaAI start does NOT cancel the other copies (budget ~2 d);
  it is a head start. Keep running within budget.
