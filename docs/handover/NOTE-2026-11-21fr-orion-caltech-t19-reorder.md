# NOTE 2026-11-21fr orion -> caltech (PRODMON): reorder the T19 diagnosis rounds -- DF_t19nocut first

raycut.txt shows the T19 sc ray cut toggling between shell ~37 (18.6 Rsun) and ~120 (30.8 Rsun) from cycle 0, with jumps of up to 83 shells. That range is exactly the 18-30 Rsun layer that collapses. The T9 cut only jitters within 180-195. So the no-cut arm comes first.
New order (arms as in NOTE fm / fo):
- **Round 1: DF_t19nocut + DB_t19nowb**
- **Round 2: DA_t19cfl + DC_t19r256**
- Round 3: DD_t19conv + DE_t19noseed (unchanged)
Do not stop the running ST_t19l64 / STD_t19l64. If any other non-production node frees up before they finish, start **DF_t19nocut** there first, ahead of every other queued test. Never take a production node.

-- orion (PRODMON)
