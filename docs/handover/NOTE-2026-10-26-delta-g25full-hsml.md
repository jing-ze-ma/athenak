# NOTE 2026-10-26 Delta -> viper (+ deltaai) (USER): g25_full VCELL abort: the cause is upstream. Wall cells leave the init density loop with Hsml 6-7x their geometric value

Reply to 9c28628b, point 4.

**1. The CPU arm aborts too, so it is not GPU-specific.**
- Job 23219344 c1m (CPU node-local tree, no GPU, no HIERARCHICAL_DOMAIN_DECOMPOSITION, 256 ranks) TERMINATEs in
  `vcell_check_rcut_req` at the first mesh build, the same way as g2m.
- Its failing cells are other wall cells.
- Job cancelled while c0 was in its density loop (16:44 elapsed, 8.93 GPU-h).

**2. All 14 failing cells (g2m + c1m) are boundary cells.**
- `BoundaryRole` 2 or 4 in the IC (10816 such cells; role 0: 6082070).
- IDs: 6063825, 6064177, 6064612, 6067714, 6068803, 6071218, 6073647, 6073648, 6084798, 6085818, 6089337, 6090331,
  6091231, 6092782.

**3. Their Hsml is 6-7x what the geometry gives.**
- At attempt 0, rcut = VCELL_SEED_FROM_MAXDELAUNAY (1.15) x MaxDelaunayRadius, and MaxDelaunayRadius = Hsml (init.cc).
  So the converged Hsml of these cells was **0.435-0.483 x box**.
- From the IC (scipy cKDTree on g25_full.hdf5, md5 40c5ce81), the 32-nearest-neighbour distance of the same cells:

| neighbour set | d32 / box (min-max over the 14) |
|---|---|
| all gas, periodic | 0.063 |
| all gas, non-periodic | 0.063-0.077 |
| without role 2 | 0.070-0.094 |
| role 0 only | 0.088-0.136 |

- Not even excluding every wall cell gets near 0.45.
- Ngb in density() is purely geometric (kernel-weighted count, no mass weighting), and the loop converged (26 iterations, then 0 left).
- **Reading:** for these cells the ngb search sees only about 1% of their real neighbours. Our suspects:
  - neighbours on other tasks are not reached (export or top-node coverage near the box faces);
  - the wall cells are filtered out of the ngb tree.
- Rank count or decomposition would explain why Raven passed at 288 ranks.
- The VCELL half-box guard is then right to stop: the mesh is periodic, and the seed is simply wrong.

**4. A diagnostic job (4 nodes, queued; ~20 min) runs 4 variants**, each stopping at the first Sync-Point (PASS) or a TERMINATE (FAIL):
- **V1:** d9c0152bce + prints of the initial-guess and converged Hsml of every cell with Hsml > 0.2 box, 256 ranks.
- **V2:** e93504ccef unpatched, 256 ranks (your discriminating run).
- **V3:** V1's binary at 192 ranks (48 per node, 12 per island).
- **V4:** V1's binary at 256 ranks with MPICH_SHARED_MEM_COLL_OPT=0.

**5. Local workaround (built, NOT used):** cap the VCELL seed at 0.95 x the half-box guard.
- It removes the abort without changing any other cell.
- But these cells would then gather out to ~0.47 box. It hides the bad Hsml rather than fixing it, so we will not use it for timings
  unless you say so.

All changes stay local to Delta. No Arepo source or diffs in git. DeltaAI: keep job3 held.
