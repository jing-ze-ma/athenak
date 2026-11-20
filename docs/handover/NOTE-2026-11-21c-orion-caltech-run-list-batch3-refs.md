# NOTE 2026-11-21c orion -> caltech (USER): batch 3 = 7 reference arms for batch 2

Follows b21808faa (batch 2). Each reference uses the same build, keys, restart, cycles and 4 ranks as its test arm; only the
feature is off. Arms: `docs/handover/caltech_runs_1121/arms_batch3.txt`. The pairs are:

| test arm (batch 2) | reference (batch 3) | difference |
|---|---|---|
| THF300_rsg | THF300_ref_rsg | `rad_m1/clean_stag_thinf=0` (it comes last on the command line, so it overrides the keys file) |
| THF100_hg | THF100_ref_hg | same |
| GCA300_rsg | GCA300_ref_rsg | build **B_v20b** = caltech-runs/v20b 694a3f9aac772158185bfee166ed77c9d7391697 (the parent of geom-ca) |
| FMP300_hg | FMM300_hg | without `rad_m1/clean_ray_flux_in=phi` |
| URF_rsg, UR8_rsg | UR0_rsg | neither extra key |
| RF300_rsg_f1 | RF300_ref_rsg | `rad_m1/clean_stag_fidp=0` |
| RF100_hg_f1 | RF100_ref_hg | same |

HGV40_12 is timing only and has no reference. Run each reference right after its test arm, or next to it, on the same kind of
node. Results go to `docs/handover/caltech_results/<TAG>/`, as before.

-- orion (for the USER)
