# NOTE 2026-10-26 Delta -> deltaai + viper (USER): diag V1-V4 at 256 ranks give no workaround; the hsmlfix test is queued but Delta's estimate is 10-28. Ask: DeltaAI tests hsmlfix-1026 without waiting for probe 3466043

Interim reply to TASK 385edc70 (hsmlfix). The full reply follows when the fix job has run.

**1. Diag job 23219733 (4 gpuA100x4 nodes, interactive, 21:14 elapsed).** Each variant stops at the first Sync-Point (PASS) or a TERMINATE (FAIL).

| variant | build | ranks | result |
|---|---|---|---|
| V1 | d9c0152bce + Hsml print | 256 | **FAIL** (VCELL `vcell_check_rcut_req`, 345 s), 29988 cells with Hsml > 0.2 box |
| V2 | e93504ccef unpatched | 256 | **FAIL** (same TERMINATE, 251 s) |
| V3 | V1's binary | 192 (48/node, 12/island) | **TIMEOUT** at the 420 s cap, still in the density loop; no verdict |
| V4 | V1's binary + `MPICH_SHARED_MEM_COLL_OPT=0` | 256 | **FAIL** (237 s), 28610 cells with Hsml > 0.2 box |

- The env setting is therefore **not** a workaround on Delta. This does not rule out 9307d33d49: that race is in Arepo's own exchange, which the env setting does not change.
- e935 does not avoid the abort either.
- The density loop converges in both V1 and V4 (26-27 iterations). So the bad Hsml is a converged value, not an iteration failure.

**2. The bad cells are a coarse lattice, nearly all on one island.**
- Converged Hsml of the bad cells: 0.20-0.46 box.
- They sit on a regular lattice of spacing box/32 (5.47e14). The IDs are the top of the ID range (>= ~6.06 M).
  - The lattice reaches the box faces and also the interior (e.g. 0.42|0.23|0.58 box), so it is wider than the role 2/4 wall cells of d6aa24a9.
  - The prints did not include BoundaryRole; `AREPO_HSML_DIAG=2` will show it.
- In both V1 and V4, every bad cell is on ranks 80-95 (10 of the 16), which is one island (16 ranks per island). No other island has any.
- **Our reading:**
  - The decomposition gives most of the sparse lattice to that island. So "one island" may reflect geometry rather than a per-island fault.
  - These cells can reach their neighbours only through exports, which depend on the top-leaf ranges of other ranks. A corrupted top-leaf exchange would produce exactly this pattern.
  - The count changes a little between runs (29988 vs 28610). That would fit a race, but it is not proof.

**3. hsmlfix-1026 on Delta: built, waiting in the queue.**
- Bundle md5 f547eaaf and `git bundle verify` OK. g2m/c1m/c0 built from 2d3967b136 with the job 23219344 recipe.
- The test job runs g2m only: 256 ranks, `AREPO_HSML_DIAG=2`, no MPICH env. The USER dropped c1m so that the job is shorter for the queue.
- It is queued on both partitions: interactive 23238294 and regular 23238295.
- **Slurm estimate: 10-28 06:37 (about 33 h).** Our fairshare fell after tonight's job, and 21 partition nodes are drained or draining.
- If it PASSes, Delta submits the 3-arm g25_full timing job (10 top steps, fixed builds) automatically.

**4. Ask to DeltaAI (USER):** please test hsmlfix-1026 now, **without waiting for probe 3466043 to reproduce**.
- Use your TASK's recipe: g2m from 2d3967b136, 2 GH200 nodes, 576 ranks then 256 if time allows, `AREPO_HSML_DIAG=2`, no `MPICH_SHARED_MEM_COLL_OPT`, 20 min or less, excluding gh082 and gh011.
- Please report:
  - PASS/FAIL per run;
  - the number of cells with Hsml > 0.2 box, and whether they are lattice cells (IDs >= ~6.06 M) on one island;
  - for those cells, the brute-force ngb count vs the tree's NumNgb (the `HSMLDIAG` summary).
- If probe 3466043 is still clean by then, a PASS will not prove that the fix works, because the bug would not be on DeltaAI. Say so in your NOTE. The fix's counts are still a useful check against the probe.
- First answer wins. If DeltaAI answers first and the fix is clearly the cure, Delta keeps its job queued for the 256-rank Delta check, unless viper or the USER says otherwise.

**SU:** diag 23219733 was 21.2 min x 16 GPUs = 5.7 GPU-h, charged x2 (interactive) = 11.3. Delta has ~974 GPU-h left.

All changes stay local. No Arepo source or diffs in git.
