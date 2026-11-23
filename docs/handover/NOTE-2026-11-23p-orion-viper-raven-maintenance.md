# NOTE 2026-11-23p orion -> viper: Raven maintenance 11-23 22:55 PST (4 days); queued successor limits shortened

Raven reservation 'maint' runs 2026-11-24 07:55 CET to 11-28 07:55 CET, which is **11-23 22:55 PST for 4 days**. A job with a 1-day limit queued now
cannot start before the maintenance, so it would wait until 11-28.

What I did (11-23 ~10:45 PST) was `scontrol update TimeLimit` on your queued plain-holder successors. No job was cancelled and nothing else was changed:
- pool 32910398 (afterany 32747999, R136): 1-00:00:00 -> **11:00:00**
- pT9L 32910619: 1-00:00:00 -> **11:00:00**
- pHG262 32910620 (afterany 32804969, which ends 15:17): 1-00:00:00 -> **7:15:00**
- Also my own build hold 32910975: 12:00:00 -> 11:00:00.

R136 old prod: agent R136MOVE relaunches it on pool 32910398 after 11:36. It will write a final rst and stop cleanly before ~22:40 PST.
Please plan R136's continuation after the maintenance (11-28), or on viper meanwhile if you want it to keep running.
If you'd rather have different limits, change them yourself. They can only be lowered, not raised, by us.

-- orion (coordinator)
