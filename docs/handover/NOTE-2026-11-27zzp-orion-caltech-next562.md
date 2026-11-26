# NOTE 2026-11-27zzp orion -> Caltech: PRIORITY 1, rerun the post-push round on caltech-next-562304fa1

Your 27zzo (thanks) tested b58f22fd0. Both of its failures are now understood:
- **ray_cut_slab** maxdiff_retry 2.39e+35: b58f22fd0 lacked the regenerated dt_retry visitors, so a retried step did not roll back the CONSERV stack's new members. d948eedaa fixed that (Orion CPU passes ray_cut_slab).
- **idort_beams** F_nocone 329: fwb_line = true lets the line-mean F error grow in vacuum when the cone cut is off. 562304fa1 sets fwb_line default false (the release behaviour); the p_ledger test names it true.

**caltech-next-562304fa1** (pushed) = d948eedaa + that one commit. Please do the 27zzo-orion list (items 1-4) on 562304fa1 instead of b58f22fd0 / d948eedaa:
1. CUDA 12.9 + g++13 build (built_in_pgens + star_vet), tiny_bsg_rel5 pass1.
2. Quick + FULL list; compare with your 27zzk (f42d79032) and 27zzo (b58f22fd0), list moved lines.
3. Combined star check REL f42d79032 vs NEXT 562304fa1, fresh 300 cycles bsg / m35mms / t9l62 (as specified in 27zzo-orion item 3).
4. The 4 production 2-step restarts on NEXT (bitwise).
If 3/4 already ran on d948eedaa, rerun only bsg (fwb_line moves it) and keep the rest. Idle test nodes only, 4 single-GPU arms per node. Post a NOTE back. Do not touch vet-clean-release.

-- orion (coordinator)
