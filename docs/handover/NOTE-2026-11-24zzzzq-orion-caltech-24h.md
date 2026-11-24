# NOTE 2026-11-24zzzzq orion -> Caltech: 24 h plan for the 4 stars: 2-node timing and projections (USER)

The USER wants all 4 stars to finish within one day; by orion's estimate one node is enough only for bsg (6 h); t9l62 (tlim 2.682e6), m35mms (8.64e5) and M35ZAMS (8.64e5) are near or over 24 h. Please measure, on idle nodes, T2 binary (7d7c4a71), settings of the U arms + opac_drag true, thinf=false, mg 4 for m35mms:
1. **t9l62, 8 GPUs (2 nodes) vs 4:** from t9l62.00015 (outburst), 100 cycles each (s/cyc over the last 50, MPI wait, dt, CW), and the same from a quiet / late rst if one exists (dt ~97): which s/cyc?
2. **m35mms, 8 GPUs vs 4:** the rst has 4 blocks (336x128x128), so 8 ranks need a fresh start with 8 blocks (e.g. meshblock nx3 = 64, layout_tile 64): say if 8 blocks work and the s/cyc of 100 fresh cycles at each layout; the dt there is 12 in both. If tile 64 changes the preconditioner so that its rise, say so.
3. **M35ZAMS, 8 vs 4 GPUs,** fresh, 100 cycles (which blocks the mesh allows); s/cyc.
4. **Projection from your own files:** hours from the T2 production histories (hst + s/cyc by cycle bin) for each star up to where it stopped, dt(t) and s/cyc(t) tables (8 bins), and t9l62's full rc3 history (dt by t to 6.7e5 and the later rst t9l62.00031): cycles and hours per bin. Post only numbers.
Priority: below the sdirk2 runs already queued (S5F done, H pair) and above nothing. One table per item.

-- orion (coordinator)
