# Handover 2026-09-22c (viper), written ~03:30

Continues HANDOVER-2026-09-22b.md.  Branch `rt-integration`, local HEAD **8e41808b** plus
uncommitted clean-up edits (see 4); fork still at **cde8940e** (nothing pushed).  Memory notes as of
now in `docs/handover/claude-memory-2026-09-22c/`; **read `master-todo-2026-09-22.md` there first**:
it holds the agreed order (finish agents -> CLEAN UP -> lists A/B/C) and the three continuation lists.

## 1. Commits since 22b (all gated; RESULTS/README in each tests_* dir)

| commit | what |
| --- | --- |
| 0eb7e2c4, 544cb36b | tabulated-EOS restarts bitwise for MHD and on the cubed sphere (`tests_mhd_rst`, `tests_cs_rst`) |
| 2268e633 | face-centred cs seam halo has NO multi-block defect (static scan; note retracted) |
| ec0a9724 | cell-centred corner slots 48-55 resampled across seams (`tests_seam_cc_corner`) |
| 3c5d4846, 2db7095c | `problem/ck_spherical`: spherical thermal two-stream for the ck path; the beam form of 3c5d4846 was WRONG and is retracted in 2db7095c (local absorption restored) |
| 13a7cf30, 916dc953 | `problem/ck_beam_sph`: pseudo-spherical direct beam (chords through shells, twilight); exact photon budget (old beam loses 14-24 %) |
| b4c6a2f2 | hipcc build fix (DualView templates in 9 files) -- the branch did not build on GPU |
| ed188f69 | **the HEAD `ck_spherical = false` path was a non-conservative HYBRID (area-weighted apply x plane-parallel sweep, from the he4 merge 1159a8f3) -- repaired to the production plane-parallel form**; `problem/ck_implicit` (backward-Euler ck column solve: gap 1e-9 but 5x cost, diverges at 10x dt; not production-ready) |
| 213095e4 | `<hydro|mhd>/rad_angular` master switch for horizontal radiative conduction |
| 131eabe8 | axial angular momentum: cs is NOT worse than sp (`tests_cs_angmom`) |
| 76359975 | per-region cs error gate (`tests_cs_regions`): vertex order 1.6; angular ppmx/wenoz cut the vertex error 3x on tests; cs_wellbalanced_src does NOT help |
| dd317c25 | switch inventory (`docs/dev/switch_inventory_2026-09-22.md`, 432 parameters, verdicts) |
| 8e41808b | GPU A/B analysis (`tests_gpu_ab_0922`) |

## 2. Findings that matter for the next production
* `ck_spherical` alone cools the night side 330-550 K at 1e-6..1e-3 bar, which pins eta at max_eta (1e13)
  in the r/Rp 1.23 shell and cuts dt 20 -> 12.5 s (Ohmic limit); `ck_beam_sph` (twilight heating)
  restores the night top to within 45 K of production.  Use the two TOGETHER.  Throughput 3.4x lower than
  the production binary; re-examine max_eta.  Isolation of the switch from the binary: arm
  `bench/ck_sph_ab/offfix` (submitted ~03:30, HEAD binary, ck_spherical = false).
* Horizontal radiative conduction is inert for dhj (noang == sph to 2e-3, not cheaper): `rad_angular = false`.
* Angular ppmx/wenoz: no vertex-localised flow exists in the real atmosphere; 1.44x cost -> not for dhj.
* Well-balanced cache: bench/cs_wb_cache every1 vs every10 (apu, 16 rot, running; death of the old WB arm
  was at rot 12.27 -- it ran WITHOUT cs_wellbalanced_src, its NOTES.md is wrong).
* `<mhd>/fofc` was never on; dhj.log's fofc column is meaningless (use eos_dfloor/efloor/tfloor/eos_fail).
* Draft next input: `inputs/production/deep_hot_jupiter_cs_prod4.athinput` (being written; DRAFT).

## 3. User decisions on the clean-up (docs/dev/switch_inventory_2026-09-22.md)
DELETE 70 abandon switches outside rad_m1 (+ bc_inflow_ghost, vdamp_bot_mean_only, mlt_relax_down_time,
wb_ramp, opac_tmin, and implicit-column modes 1/2 with rt_impl_tau_min/dtmax/tau_blend); KEEP
rt_layer_legacy; rad_m1 untouched (paused).  DEFAULT ON: ck_spherical and ck_beam_sph (both refused on
Cartesian), rt_use_cons, rt_col3_ex_iter (check inert on RG/dhj), mlt_split_deposit, rt_split.
wb_cache_every stays 0.  Not decided: the 36 UNSURE (kept).

## 4. In flight at writing
* Sonnet clean-up agent: groups G1 (conduction), G2 (two-stream/ck/dhj), G4 (red giant), G5 (WB/cs/polar):
  387 lines removed so far, uncommitted, gates in `tests_cleanup_0922/`.  G3 (box_convection, 37
  switches) and G6 (column modes 1/2) wait for Opus (two 529 failures).  Then the default flips.
* GPU: bench/cs_wb_cache every1/every10 (rot 1.4 at 03:00, no collapse); bench/ck_sph_ab/offfix.
* Sonnet: prod4 draft input + audit of bench NOTES.md vs inputs (`tests_bench_audit_0922/`).

## 5. Standing
M1 paused (`m1-paused-back-to-dhj`).  Delegation: Sonnet for mechanical/gated work, Opus for judgement.
The `origin` remote URL embeds a token (user action).  cs_mhd_prod3/rst = 46 GB (retention rule: user).
