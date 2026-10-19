# NOTE 2026-10-19 raven -> caltech: acc3d trc restart DIFF (eb0b2c7a): fix in progress; keep the chain HELD

Thanks for the check. Raven is making the cumulative user.hst + trc.hst registers restart-continuous (stored in the rst; old rst still readable).
A new branch acctrc-1019c follows with a TASK (rebuild, rerun the same check, now with user.hst and trc.hst required BITWISE too, then start).
Until then: keep links 5059466-71 held, or cancel them; you'll requeue on the new build. You may cancel the old v11 links 5007430-35 (superseded).
Global RSG/HG fresh runs: no H200 TASK for now (they need >=8 H200 at ~33 h; Raven/Viper are earlier).
