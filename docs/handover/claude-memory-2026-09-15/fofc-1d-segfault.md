---
name: fofc-1d-segfault
description: "Pre-existing crash on wb_allmach_cubedsphere — FOFC + 1D hydro + PLM segfaults, unrelated to the general EOS work"
metadata: 
  node_type: memory
  type: project
  originSessionId: d0224e66-3587-4359-a545-9608682177c3
  modified: 2026-08-07T16:00:39.420Z
---

On branch `wb_allmach_cubedsphere` (and inherited by `general-eos`), running the 1-D Sod
shock tube with `<hydro>/fofc = true`, `reconstruct = plm` and `nghost = 3` **segfaults**
(exit 139, signal 11), for both `eos = ideal` and `eos = general`.

Verified pre-existing on 2026-08-07: stashed all general-EOS changes, rebuilt at HEAD, ran
the unmodified ideal-gas binary — same segfault. It is **not** caused by the general-EOS
work.

2-D MHD Orszag-Tang with FOFC runs fine on both EOS paths, so it looks specific to the 1-D
hydro FOFC path. Not yet debugged. Note `hydro_fofc.cpp` guards its j/k blocks with
`multi_d`/`three_d`, so a missing dimensional guard is not the explanation.

Related: [[general-eos-project]].
