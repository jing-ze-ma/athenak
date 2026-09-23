---
name: overnight-rules-0915
description: User asleep from 09-15 ~11:00 - NO watches, NO reactions on the production chains (prod_w4 11712113-5, box_w4 11712116-8): do not cancel, restart, resubmit or edit them; only record the results of the three agents still running (hybrid PCR port, bottom inflow, handover addendum) in memory; the user checks the chains themselves
metadata:
  type: feedback
---

The user said: "don't monitor anything because then you have to react. I'm going to sleep. can I trust you with this?" Both startup watches were stopped. Standing behaviour until the user returns: no Monitor/background waits; if an agent reports, write memory only; no new jobs beyond the agents' own small gate arms; never touch the production directories or jobs; when the user returns, give the state from squeue + out.txt in one check. See [[no-heavy-monitoring]], [[save-tokens-everywhere]].
