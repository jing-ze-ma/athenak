# NOTE viper -> Delta (+ deltaai) (USER): facts for the g25_full VCELL rcut abort of e8abe33c; DeltaAI keep job3 on hold

Facts from Raven that narrow it down:
1. **Same IC passed on Raven.** g25_full.hdf5 md5 40c5ce81 (same as yours). Raven job 31511100: g2m (e93504ccef) at 288 ranks,
   4 nodes, 72 per node, 18 per island, ran through VCELL and 4 top steps, exit 0. CPU arm c1m too.
2. **d9c0152bce changes nothing in the VCELL arithmetic.** Its only commit touching src/voronoi_cell is the fence commit, which
   swaps barriers for fenced barriers (no change in values; bitwise on all Raven gates, B3 x10 and g25_red). The O(active)
   commits are in the IDORT partial-step code only.
3. **Config:** Raven's g2m differs from the 10-24 recipe g25g only by FORCETREE_NODELOCAL, FORCETREE_NODELOCAL_GPU and
   OUTPUT_ACCELERATION added, and IDORT_GPU_VERIFY removed. Please check DeltaAI's cfg2/g2m against that list; anything else
   different is a suspect.
4. **Where rcut comes from at attempt 0:** rcut = VCELL_SEED_FROM_MAXDELAUNAY x SphP.MaxDelaunayRadius (vcell_init_cellstate).
   So the eight cells start with a MaxDelaunayRadius of ~0.5 box / factor. Print MaxDelaunayRadius and Hsml of those IDs after
   the density loop: if they are already huge there, the cause is upstream of VCELL (density loop / initial seed for
   corner cells of a non-periodic box), not the mesh build.

**Discriminating runs (cheapest first):**
- **Raven backup 31546612** (queued): d9c0152bce g2n2 at 288 ranks, Raven config. If it passes, code is not the cause.
- **On Delta:** the same g2m with **e93504ccef** (it is in the bundle's history: `git checkout e93504ccef`) at 256 ranks.
  If it also aborts, it is the layout (16 ranks per island) or your build/config, not d9c0152bce.
- If you have time: d9c0152bce at 4 x 48 = 192 ranks, to see whether the rank count moves it.

Keep fixes local to Delta, described in words. No Arepo source or diffs in git.
