# NOTE deltaai -> viper + Delta (USER): ack TASK 385edc70; probe resubmitted WITHOUT MPICH_SHARED_MEM_COLL_OPT=0 as 3466302

- Ack TASK-2026-10-27-viper-deltaai-hsmlfix. The bundle arepo_hsmlfix_1026 is **not on DeltaAI yet** (waiting for the user to copy it).
- **The problem with 3466043:** jobh.sub exported `MPICH_SHARED_MEM_COLL_OPT=0`. That is the DeltaAI workaround for the same race class that 9307d33d49 fixes.
  - With it on, the probe could not reproduce the bug.
  - So a clean result would NOT have meant "the bug is not on DeltaAI".
- **USER go:** cancelled 3466043 (never ran, 0 SU) and resubmitted as **3466302**.
  - Identical to 3466043 except `MPICH_SHARED_MEM_COLL_OPT` is unset; the job log prints its value to prove it.
  - Same binary (d9c0 + HSMLPROBE print, md5 31b1e582).
  - Same layout: 2 GH200 nodes; n576, then n256; 12 min; ~3 SU.
  - Excludes gh082 and gh011.
- **Risk:** without the setting it may hang (as the island sort did) rather than print. The per-run timeouts cap it, and we'd report a hang as "reproduces something".
- **Next:** we follow your decision rule on 3466302's result. If it reproduces, the fix test also runs with the setting unset, as you asked.
- job3 stays cancelled.
