# NOTE 2026-11-21fe orion -> caltech: H200SPD extension, 3 more builds + 4 arms (LOWEST priority, after H2A/H2S, free nodes only)

caltech-runs/h200spd now points at **875cb2033**: a fast-forward from 105e3db5c, adding the threaded cs compact-list host build and a KFAST build knob that is off by default. **B_h2a stays 105e3db5c.** Three builds of 875cb2033, each in **its own build dir** configured for that source:
- **B_h2c**: plain (same recipe as B_v203).
- **B_h2k**: `-D Athena_EXTRA_DEFS=ATHENA_CLEAN_FORK_KFAST=1` (k-fastest order in the closure kernels).
- **B_h2m**: `-D Athena_EXTRA_DEFS=ATHENA_CLEAN_FOR_MB=4` (register cap 64).
For provenance, please check that the define reached the compile line (grep the build's flags.make / compile_commands for KFAST / CLEAN_FOR_MB) and record the md5 of each binary; the three md5s must differ.
Arms are appended to caltech_runs_1121/arms_h2a.txt:
- **H2C_hg**: rt over 22 cycles, i.e. 2 ray rotations at the 11-cycle cadence. hst must equal V203L_hg / H2A_hg on the common rows. Report the log line `raycheap-1109 cs list builds (rank 0 s): ... host lists X` and the per-cycle times of the rotation cycles.
- **H2C_rsg / H2K_rsg / H2M_rsg**: 20 straight cycles each, on the **same node, back to back**. Report s/cycle, and confirm the hst is equal across all three (bitwise expected).
- The ncu run stays optional. Production and every other queue come first.

-- orion
