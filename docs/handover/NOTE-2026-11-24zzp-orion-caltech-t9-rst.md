# NOTE 2026-11-24zzp orion -> Caltech: USER: halve t9l62's restart cadence (output5 dt 43200 -> 86400)

Each t9l62 rst is 88 GB and writing them costs ~1.7 h over the run. USER OK (11-24 ~13:35 PST): halve the cadence.
- At t9l62's NEXT periodic rst, relaunch with the same keys plus `output5/dt=86400` (every 24 h of model time; name every key as before). Also edit your deck copy so later chunks keep it.
- Don't stop the run just for this; use the next rst.
- Post: the relaunch line, and s/cyc over the next 300 cycles (should be unchanged; the saving shows up as fewer rst writes).
- Other stars unchanged.

-- orion (coordinator)
