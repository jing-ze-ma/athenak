---
name: gpfs-quota-wall
description: /orion/u (= /u/jinma) has a hard per-user GPFS quota that silently TRUNCATES simulation output instead of aborting the run
metadata: 
  node_type: memory
  type: project
  originSessionId: b25dc958-761c-420e-9c05-1aabb20efebb
  modified: 2026-08-23T00:00:00.000Z
---

`/orion/u` and `/u/jinma` are the SAME gpfs filesystem (`orion_u`, 14 T, ~11 T free). `df`
therefore tells you NOTHING about whether you can write — there is a hard **per-user
quota**, and it is the binding constraint. Hit on 2026-08-14 at ~113 GB of usage.

## The dangerous part: failures are SILENT

- Writes fail with `OSError: [Errno 122] Disk quota exceeded` in Python, but **AthenaK does
  not abort** — a running job keeps computing and logs no error at all
  (`grep -icE "error|quota|fail" log.out` returned 0). It just writes SHORT `.bin` files.
- Three dumps were corrupted before this was noticed: two truncated (68 MB and 34.6 MB of
  an expected ~70.8 MB) and one that was **216 bytes**.
- **Detection method:** AthenaK `.bin` dumps from one run are all the same size to within a
  few TENS of bytes (the header embeds the time as text, so it wobbles: 70796050 /
  70796051 / 70796065 for the same run). A truncated dump is unmistakable — orders of
  magnitude off and landing on a round binary boundary (68157440 = 65 MiB,
  34603008 = 33 MiB, or 216 bytes for a dump that got nothing).
  `ls -l <run>/bin | awk 'NR>1{print $5}' | sort | uniq -c | sort -rn` shows the
  distribution at a glance. Do this before trusting any plot or comparison built on dumps.
- `dd if=/dev/zero of=... bs=1M count=20` is the quick write test: it reports
  `0 bytes copied` when over quota.

## Reading the quota

`mmlsquota` is NOT installed on the login nodes, and `/usr/lpp/mmfs/bin/mmlsquota` fails
with "GPFS is down on this node". `quota`/`lfs quota` do not apply. So the limit cannot be
read from where I work — infer it from behaviour, or ask the user to run `mmlsquota -u
jinma orion_u` somewhere that works.

## Budgeting

Usage as of 2026-08-14: `ATHENAK` 49 G (of which `run/` ~40 G), `AREPO` 32 G, `SLH` 26 G,
misc 2.2 G. A 96 x 192 x 192 `hydro_w` dump is **~71 MB**, so a 21-dump run costs ~1.5 GB
and a PAIR of them ~3 GB. Budget this BEFORE submitting high-resolution runs — check
headroom with a 200 MB `dd` test, not with `df`.

Freed 11 GB by deleting `run/{powerlaw_test,radiative_test}/bin` at the user's suggestion;
see [[run-directory-untouchable]] for why only `bin/` and not the whole directories.
Related: [[solar-convection-general-eos]], [[scratchpad-not-visible-to-compute-nodes]].

## Hit again 2026-08-23, and the cheapest thing to reclaim

`git fetch` died with `fatal: write error: Disk quota exceeded`; the `dd` probe wrote
0 bytes. Usage was ~108 GB of top-level home (ATHENAK 48 G, AREPO 32 G, SLH 26 G,
conda-envs 1.8 G) — so the wall sits near ~110 GB, consistent with the 2026-08-14 hit.

The reclaim that cost nothing real: **core dumps**. Four of them sat in
`ATHENAK/eostest/{head_sod,t_out,amrbase_out,diag_out}/core.orion02.*` — 57 GB
*apparent*, ~1.2 GB of actual blocks (they are sparse, so `du -sh` on the parent hides
them; `find ... -name 'core.*' -printf '%s %p\n'` is what surfaces them). Deleting them
freed enough to fetch immediately. GPFS quota counts real blocks, not apparent size.

So when over quota: `find /orion/u/jinma -maxdepth 4 -type f -name 'core.*' -size +50M`
BEFORE proposing to delete anything the user actually cares about. `run/` is off limits
([[run-directory-untouchable]]) and the sibling data dirs (eostest, soleos, idflr) are
simulation output, not junk.
