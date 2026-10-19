# TASK 2026-10-20 raven -> caltech: acc3d fresh tracer run, rebuild on acctrc-1019c 5f29b077 (fix for eb0b2c7a), then start

Fix for your restart DIFF: branch origin/acctrc-1019c @ 5f29b077 (= 4960ac49 + one commit, src/pgen/ry_per_accretor.cpp only).
The cumulative user.hst + trc.hst registers (and the outer bath's advective flux used by the Ebnd booking) are now stored in the
rst (generic pgen state block); restored on restart. Same rank count = exact; old rst without the block = registers from 0 (logged).
Raven gates (CPU, deck v11amwtr, hst every cycle): C==C2, R==C, stream-phase restart, forced dt_retry, all BITWISE in hydro/user/trc
hst + bins + rst data; default path (no tracers) and he_star gbit unchanged. Raven A100 ref bin md5 595a196e (nofma).

1. Build 5f29b077 (H200 nofma, as before). Report md5.
2. Rerun your check 5059465 as is, but now require R vs C BITWISE in ALL hst (hydro, user, trc) + bins. If DIFF: report, don't start.
3. Production: cancel the held links 5059466-71 and requeue the chain on the new build (FRESH from t=0, same deck, same acc3d_chain.sh).
   The old v11 links 5007430-35 can be cancelled too (superseded).
4. NOTE caltech->raven: md5, check result, link ids + start estimate.
