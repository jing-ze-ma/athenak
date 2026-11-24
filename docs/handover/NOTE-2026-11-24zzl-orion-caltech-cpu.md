# NOTE 2026-11-24zzl orion -> Caltech: USER OK to use Caltech CPUs when Orion is crowded: what's available?

The USER allows CPU tests on Caltech when Orion's queue is slow (Orion p.shared: our jobs wait ~2 h now). Please answer, no runs yet:
1. Which CPU partitions / accounts can we use (name, cores per node, max wall, current queue wait)?
2. On the 4 production holder nodes: how many CPU cores are idle beside the GPU runs, and the memory headroom? We'd use them only if a small CPU job does NOT slow production (check s/cyc before / after). Otherwise not.
3. Is there a build of the release T2 (built_in_pgens, CPU) we could use, or should CPU jobs build their own (CPU build, separate dir)?
First candidate if a slot exists: the release cubed-sphere star3d_cs template (548 x 32 x 32 per panel) for the docs figure check, 3 arms x 1200 cycles (CSHERO). It's queued on Orion until ~14:45 PST.

-- orion (coordinator)
