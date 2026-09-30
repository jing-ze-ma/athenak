---
name: test-output-location
description: All test/analysis output goes in /orion/ptmp/jinma/Athenak/ -- a separate GPFS with petabytes free, NOT the quota-bound /orion/u
metadata:
  type: feedback
---

**Save all test-related files -- run output, dumps, analysis products -- under
`/orion/ptmp/jinma/Athenak/`.** Instructed by the user 2026-08-23.

**Why:** `/orion/ptmp` is a DIFFERENT filesystem from `/orion/u` (2.8 PB, ~1.1 PB free)
and is not subject to the per-user quota that keeps stopping work on `/orion/u`
([[gpfs-quota-wall]]). It is also real GPFS, so MPI-IO works -- unlike `/tmp`, where
ROMIO cannot open AthenaK's output files at all ([[exo-fms-ck-tables]]).

**How to apply:** `-d /orion/ptmp/jinma/Athenak/<runname>`. The directory already exists
and holds earlier work (`cart_blast`, `pole_blast_*`, `sph_blast*`). Do NOT stage runs on
`/orion/u/...` "just for a quick test", and do not use `/freya/ptmp` for AthenaK work on
orion -- that path is the user's 25 TB Arepo area and is not the right home for this.
Source, inputs and the data tables still live in the repo on `/orion/u`; only OUTPUT moves.
