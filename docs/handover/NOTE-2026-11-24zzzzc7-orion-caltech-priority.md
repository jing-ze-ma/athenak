# NOTE 2026-11-24zzzzc7 orion -> Caltech: priority and trimmed queue (18:05)

Production (bsg, m35mms, M35ZAMS on 95813cf85; t9l62 on c4) first, always. Then, one item per idle node, top item first:
1. **sdorder5 (zzzzb5):** CUDA build + GPU suite (full list) + old-rst check. Holds up the merged release.
2. **cs pattern (zzzzc6 section 1):** Bp, E, G, H (4-5 min each, together ~20 min on one node).
3. **PRF_zam** (zzzzc2 section 1, M35ZAMS): decides whether the ring-sum fix of SPEED24 pays there. One arm.
4. CPU-only, any idle core, never on a production node: vceil cells (zzzzc4), recurrence (zzzzc5).
**Dropped (don't run):** PRF_t9o, PRF_m35c (the two profiles we have are enough), MLQ0 / MLQ1 / MLQ1r and section 3 (LW) of zzzzc2 (merging gains ~1.1 h only in the quiet phase; the running t9l62 will show the outburst and late cost for free), any repeat of HG arms. If an arm of these has started, let it end and post it.
**Coming:** one NOTE with a SPEED24 A/B binary (atomics + ring sum + sum-rule dt + line-solve key) once its Orion gates pass: keep one node free for it when this queue drains; don't start anything speculative on it.

-- orion (coordinator)
