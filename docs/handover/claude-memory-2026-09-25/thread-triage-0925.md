---
name: thread-triage-0925
description: User 09-25 ("we opened too many things"): keep only the dhj production path; side threads finish their current step and stop; mg phase 2 parked; ck-fast2 waits
metadata:
  type: feedback
---
Applied 09-25 ~04:30. KEEP (production path): ck-nq2 (GPU jobs -> B vs blend -> merge), ck-hitemp (merge after), dhj-cksph-guard
(merge), ck RCE initial profile, radial A/B a320; SPARC sponge runs start automatically once those binaries/IC are in.
FINISH THEN STOP: SC-on-sp prototype (report), m1-fast3 (gates -> merge). m1-mgfuse: phase 1 (kernel fusion) only;
phase 2 (global coarse correction) PARKED. ck-fast2 (nq2 speed + scaling) NOT started until the production binary exists.

**Why:** the user felt too many threads were open.
**How to apply:** don't open new side threads without asking; propose, and prefer finishing/merging over starting.
Related: [[sparc-sponge-campaign-0925]], [[ck-nquad2-production]], [[m1-precond-mg]].
