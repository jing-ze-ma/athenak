---
name: use-rt-integration-branch
description: "STANDING (user 09-15 evening): after the merge, rt-integration (bench/wt_merge) is THE working branch for all box / red-giant / hot-Jupiter work; rg-box-pcr, polar-average-perf, lhlld, implicit-transverse-raddiff are merged into it and frozen"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: e67ee046-aa0b-46f4-881b-84099bfda99a
  modified: 2026-09-20T22:56:38.909Z
---

**Why:** the user wants one clean branch carrying every update from the branches used for rg and dhj; the wholesale header copy showed how divergent hand-ports go wrong.
**How to apply (09-16 09:20):** rt-integration b775a270 carries ALL work (labs merged, default-off switches; HIP configure needs the openmpi_gpu bin on PATH and -D CMAKE_PREFIX_PATH=/mpcdf/soft/RHEL_9/packages/x86_64/rocm/6.3.4); new commits, builds, binaries and worktrees start from rt-integration (pushed to fork). Do not commit on rg-box-pcr / polar-average-perf any more; if a production still runs on an rg-box-pcr binary, switch its next link to an rt-integration binary only after a bitwise gate. See [[rt-integration-branch]].

09-17 19:00 (user decision): the MAIN CHECKOUT /viper/u2/jinma/ATHENAK/athenak is now on
rt-integration (728e60d5 = rt-integration + the two 09-17 handover-doc commits merged in);
implicit-transverse-raddiff DELETED (it had no code beyond rt-integration); the wt_merge worktree
(with its 26 build dirs) REMOVED. Production binaries live as athena_pinned copies in the run dirs.
Global He4 work stays on he4-presn-global in bench/wt_he4 (parent rt-integration).

MERGE GATES 09-21 (tests_gate_merge/, untracked; e688efd4 vs 8841786d): box G1 (box_w8 16x16,
50 cycles, serial CPU, modes 3 and 0) BITWISE on all 5 files; red_giant_cs 100 cycles GPU is
NOT bitwise (cycle 0 and column.txt identical; dt 3.9e-4, mass/E/KE equal to 6 digits) = the
spherical two-stream now acts on every red_giant run (rt_plane_parallel is a code flag set only
by box_convection, not an input). Old RG runs (prod11) are therefore not reproducible bit for bit
with post-merge binaries. New GPU red_giant binary: build_gpu_rg/src/athena.

09-21 (user decision, SUPERSEDES the line above): he4-presn-global (12d3644f) merged into
rt-integration as e688efd4 (clean, src identical to the He4 branch; PUSHED to fork 09-21). From now on ALL
work, He4 included, is done on rt-integration in the main checkout; he4-presn-global is frozen.
bench/wt_he4 still held uncommitted edits (two_stream_rt.hpp, red_giant.cpp, he4_presn_cs/sp
inputs = the rt_kappa_hsmooth work in progress) at merge time. User 09-21: do NOT use those
uncommitted in-progress edits; rt-integration stays at the committed state (they remain only in
wt_he4, untouched).
