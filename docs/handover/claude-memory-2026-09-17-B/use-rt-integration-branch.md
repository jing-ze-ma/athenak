---
name: use-rt-integration-branch
description: STANDING (user 09-15 evening): after the merge, rt-integration (bench/wt_merge) is THE working branch for all box / red-giant / hot-Jupiter work; rg-box-pcr, polar-average-perf, lhlld, implicit-transverse-raddiff are merged into it and frozen
metadata:
  type: feedback
---
**Why:** the user wants one clean branch carrying every update from the branches used for rg and dhj; the wholesale header copy showed how divergent hand-ports go wrong.
**How to apply (09-16 09:20):** rt-integration b775a270 carries ALL work (labs merged, default-off switches; HIP configure needs the openmpi_gpu bin on PATH and -D CMAKE_PREFIX_PATH=/mpcdf/soft/RHEL_9/packages/x86_64/rocm/6.3.4); new commits, builds, binaries and worktrees start from rt-integration (pushed to fork). Do not commit on rg-box-pcr / polar-average-perf any more; if a production still runs on an rg-box-pcr binary, switch its next link to an rt-integration binary only after a bitwise gate. See [[rt-integration-branch]].
