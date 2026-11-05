# NOTE 2026-11-05 Caltech -> viper: Raven maintenance (USER)

USER: "raven is going to be maintained soon".
1. Caltech KEEPS its R136a1 128² runner armed (it starts on a free Caltech node, ~11-06 04:00). The Raven R136a1 probe is
   still welcome if it finishes before the maintenance, but it is not the primary run.
2. Please report the exact Raven maintenance window (start/end, which partitions/filesystems, /raven/ptmp access) and
   whether the RSG E / HG gap plan (11-07 ~05:00 .. ~11-10/11) is still feasible. If the window overlaps, say so,
   so the user can decide before copying ~140 GB of restart files.
3. Do not start any chain that would be killed by the maintenance without a final rst; plan links so the last one ends
   at a rst before the window (touch FINAL / shortened TimeLimit as in your reverse plan).
