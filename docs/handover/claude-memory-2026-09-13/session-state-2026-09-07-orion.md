---
name: session-state-2026-09-07-orion
description: Entry point after the viper->orion handover of 2026-09-07; what was pulled, what is pending, what to do first
metadata:
  type: project
---

On 2026-09-07 the user updated `polar-average-perf` on the fork from viper and asked me to
pull it here on orion. Pulled by fast-forward to 5e9b1273 (stash/pop of orion-local edits,
which were COMMITTED AND PUSHED at the end of this session (db086512, 60c278c7, 802bf827): data/exo_fms_ck/PROVENANCE.md,
both dhj inputs (rst output block + meshblock timing comments), src/pgen/solar_convection.cpp).
The orion-local CLAUDE.md was replaced by the tracked one; the old copy is gone with the scratchpad.
The ck RT tables (data/exo_fms_ck/{ck,cia,CE_tables}) are gitignored and untouched.

Handover documents: docs/handover/HANDOVER-2026-09-07.md (read first), docs/handover/scripts/,
docs/handover/claude-memory-2026-09-07/ (imported into this memory dir, see the second half of MEMORY.md).
Viper is in maintenance 2026-09-07 12:00 -> 2026-09-12 12:00; its running jobs' results are only readable from viper.

**Why:** the two machines cannot see each other's run dirs; the handover doc is the only shared state.
**How to apply:** next steps in order: (1) CPU build with PROBLEM=deep_hot_jupiter_rt, (2) run
tst/test_suite/rad/test_rad_dhj_ck_cpu.py and test_rad_cs_raddiff_cpu.py, (3) any GPU build must
include d3d74f2b, (4) cs MHD determinism test (restart twice from one rst, 200 cycles, bitwise diff),
(5) regenerate eos_table.txt via <mhd>/eos_table_dump before any T/p plot. See [[cs-mhd-prod-nan-rot41]],
[[wb-hydrostatic-scheme-cs]], [[radiative-conduction-deep-interior]]. Orion memory of the same period
lives in [[session-state-2026-09-04]] and [[test-output-location]].

**Update, 2026-09-07 evening:** steps 1-2 DONE. Both tests PASS on orion at 5e9b1273
(test_rad_dhj_ck_cpu 45 s, test_rad_cs_raddiff_cpu 36 s). Persistent CPU build (serial+OpenMP, no MPI)
is at build_dhj_cpu/src/athena. Toolchain that works on the orion login node:
`module purge; module load gcc/13 cmake/3.28 anaconda/3/2023.03; export CC=gcc CXX=g++`.
TRAP: system python3 is 3.6 and run_test_suite.py dies with `Popen ... unexpected keyword 'text'`;
use anaconda/3/2023.03. Piping `module load ... | tail` runs it in a subshell and loads nothing.
Running `python3 -m pytest <file>` from tst/ skips the driver's extra full build; tests build their own.
Remaining: steps 3-5 (GPU build with the NaN guard, cs MHD determinism test, eos_table_dump).


**End of session: pushed to the fork as 8df604da.** Four commits on polar-average-perf:
db086512 the exo_fms_ck upstream-defect audit, 60c278c7 the two dhj inputs (measured CPU
decomposition, the <output3> rst block, ck keys on ideal_xe), 802bf827 the solar_convection
restart-gravity fix (identical content to ef9561e2 on fix/solar-convection-restart-gravity,
which is NOT merged here), 8df604da docs/handover/RESULTS-2026-09-07-orion.md answering the
viper handover. The working tree is clean apart from the untracked run/ directory.
