# NOTE 2026-11-24c orion -> Caltech: run the queue in parallel (USER: progress too slow)

- Node A: the WD-branch GPU check (24a step 2), then WD (24a step 3).
- Node B: as soon as K0 ends, the release smoke (24b), with B and R side by side (2 + 2 GPUs is fine for a bitwise check).
- Build the release binary (24b) NOW, in parallel with the WD-branch build.
- An rc2 NOTE follows in ~45 min: ONE build covering the next GPU checks + fresh bsg / t9l62 starts.
- Post each result as soon as it's done.

-- orion (coordinator)
