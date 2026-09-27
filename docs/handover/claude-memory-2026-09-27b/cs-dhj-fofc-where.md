---
name: cs-dhj-fofc-where
description: 2026-09-13 per-cell FOFC map on the cs hot-Jupiter MHD (bench/cs_mhd_prod3_fofc/deathmap/analysis): FOFC fires ONLY at the top of the atmosphere (98% at i>=98, 63% in the top 10 rows, 0 below i 85), on the cold NIGHT-side hemisphere (12 of 24 blocks; 100% of flagged top columns at dfloor), NEVER on low-beta faces (0 of 1.7e9 flags with beta<1; the beta<0.5 cells are 7 per dump at i 47-52, the H2 front) => FOFC cannot substitute for cs_lowbeta_fallback; they act on disjoint regions. 3.7% of cell-stage updates redone at 1st order for +4.8% dt; efloor 8.4x cs_mhd_prod3
metadata:
  type: project
---

Instrument: flag totals reconcile with dhj.log fofc to 0.13% (ghost_zones=true dumps). Seam/vertex enrichment only
1.1-2.5x (85% of flags off-seam). dfloor still exceeds fofc 1.7x: FOFC keeps rescuing the same vacuum slab at the
top rather than removing it. The cs_mhd_prod3_fofc arm (fallback off) tracked cs_mhd_prod3's energetics and dt but
died at rot 11.4 (unexplained, non-reproducible). **How to apply:** keep cs_lowbeta_fallback ON for cs dhj; FOFC is
optional and only touches the night-side top; the real fix for the top slab is the density/energy floor treatment
there (dfloor 5e-14 cells at i>=118), not FOFC. See [[cs-wb-arm-and-ke-gap]], [[fofc-gpu-verified]].
