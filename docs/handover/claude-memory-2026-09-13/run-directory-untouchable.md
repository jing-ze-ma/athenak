---
name: run-directory-untouchable
description: "The run/ directory holds ~40GB of simulation data and must never be committed, cleaned, or modified"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: d0224e66-3587-4359-a545-9608682177c3
  modified: 2026-08-14T16:12:10.634Z
---

`run/` in the AthenaK repo is ~40 GB of untracked simulation output and analysis scripts. Never `git add` it, never `git clean` it, never modify its contents. When doing any git operation that could touch the working tree (checkout -f, reset --hard, clean), verify `run/` survives afterward.

**Why:** It is irreplaceable simulation data that is not backed up in git, and a stray `git add -A` at the repo root or a `git clean -fd` would either bloat the repo or destroy the data.

**How to apply:** Stage with explicit paths (`git add src/ ...`), never `git add -A` at the repo root. Scope any `git clean` to a subdirectory. `git checkout -f` and `git reset --hard` are safe for `run/` since it is untracked, but confirm with `du -sh run` afterward. See [[athenak-fork-git-setup]].

**The rule is about ME not touching it on my own initiative — the user can and does grant
exceptions.** On 2026-08-14, at the user's explicit suggestion, `run/powerlaw_test/bin` and
`run/radiative_test/bin` were deleted (5.5 GB each, 11 GB total) to clear a GPFS quota wall
(see [[gpfs-quota-wall]]). Everything non-bulk was KEPT in both dirs — athinput, the pinned
binary, submit script, logs, and crucially the complete `.hst` history files (101 rows to
t=20000), which hold the mass/energy/KE time series the writeup rests on. So both control
runs remain fully reproducible and their quantitative results survive in
[[solar-convection-test]]. When asked to free space, prefer deleting only `bin/` for this
reason rather than whole run directories.

## Exception: user-directed cleanup (2026-08-23)

The rule is "do not clean it on your own initiative", not "never touch it". When the user
explicitly asks for a cleanup, the safe categories — and the ones they approved — are:

- **`.athdf` files are caches, not data.** `plotsun.py` and `plotcool.py` both call
  `bin_convert.write_athdf(...)` unconditionally right before reading, so every `.athdf`
  is regenerated on the next plot. Deleting them is free. 201 of them in
  `run/sun_test` + `run/cooling_test` were 5.1 GB.
- **Superseded `.rst` files.** Keeping only the newest restart per run directory freed
  1.34 GB across `ATHENAK/soleos` + `ATHENAK/eostest` (1373 of 1411 files).
- **Core dumps** — see [[gpfs-quota-wall]].

What the user declined to delete, so do not offer it again unprompted: the `prev_179*`
parameter-scan runs in `run/sun_test` (their `.bin` dumps are not regenerable), the
`sun_test` core dump, and either half of `eostest/refactor_chk` vs `refactor_after`.
Those two trees are 47/50 byte-identical by sample, but `amr_general`, `smr_general` and
`smrm_general` DIFFER — that difference is the refactor comparison itself.

Always propose the list and get an explicit choice before deleting; never infer scope.
