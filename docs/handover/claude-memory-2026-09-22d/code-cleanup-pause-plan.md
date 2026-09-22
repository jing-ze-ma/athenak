---
name: code-cleanup-pause-plan
description: User 2026-09-22 - after the running agents finish, pause feature work and clean up the code; many switches should be abandoned, some made default-on
metadata:
  type: project
---

User, 2026-09-22 ~02:00: once the running agents finish (implicit ck column solve, cs per-region gate)
take a PAUSE from new features to clean up: "some of the switches should just be abandoned and some
should be default on".

**Why:** the branch has accumulated hundreds of default-off switches (rt_*, rad_*, ck_*, wb_*, implicit_*,
dbg_*); tonight a merge regression hid behind one ("ck_spherical = false" at HEAD is a non-conservative
hybrid nobody intended), agents' bitwise-default discipline keeps dead paths alive, and failed
experiments were kept "as options" (implicit_accel = anderson, implicit_trans_limit, dbg_trans_memory,
dbg_tensor, implicit_closure_relax, ...).

**How to apply:** do NOT start deleting on my own. First produce an INVENTORY (read-only agent): every
runtime parameter of src/utils/two_stream*.hpp, correlated_k, src/diffusion/conduction*, src/rad_m1,
the dhj / box_convection / red_giant pgens, WB (hydro/mhd), with default, where read (file:line), what
production inputs set it to (bench/cs_mhd_prod3, sp_mhd_prod3, bstar_fecz prod_w7, hestar_fecz box_w8,
RG prod11), which gate/test covers it, and a proposed verdict: ABANDON (failed experiment / superseded /
never used), DEFAULT ON (the fix is the correct physics: e.g. ck_spherical, wb_cache_every = 1,
implicit_gas_newton, f_source = wb, rad_angular = false for dhj?), KEEP AS OPTION. The user decides per
switch. Changing a default breaks bitwise gates by design: regenerate the reference numbers in the same
commit and say so. Candidates I already believe: default-on ck_spherical (+ fix the off path to the
production plane-parallel form or remove it), ck_beam_sph, MHD/cs restart caches (already unconditional);
abandon the HEAD hybrid, implicit_accel, dbg_trans_memory, implicit_closure_relax(_thin). Related:
[[dhj-improvement-list-0921]], [[m1-paused-back-to-dhj]].

## User decisions 09-22 ~04:15
- DELETE the 70 abandon switches outside rad_m1 (agent running, tests_cleanup_0922/); KEEP rt_layer_legacy;
  rad_m1's 31 abandon switches stay for now (M1 paused).
- DEFAULT ON (queued after the deletion, same files): ck_spherical (true on cs/sp, REFUSED on Cartesian),
  ck_beam_sph (true, REFUSED on Cartesian), rt_use_cons, rt_col3_ex_iter (check bitwise-inert on red giant
  and dhj), mlt_split_deposit, rt_split (inert; delete the off path).
- wb_cache_every stays 0 (exact); forbid/warn N > 1. rad_angular = false goes into the dhj INPUT only.
- Not decided: the 36 UNSURE (list printed 04:05) and the M1 default-ons.
- GPU A/B verdict (tests_gpu_ab_0922): keep ck_spherical + ck_beam_sph TOGETHER; sph alone cools the night
  top 500 K and hits the Ohmic dt cap (max_eta 1e13) -> 2.8-3.4x throughput cost, to be re-examined via
  max_eta; noang inert -> dhj input; ppmx/wenoz useless for dhj.
DONE 09-22 06:00: clean-up committed (47eb060e, 912ef43c, ac42334d; -4779 src lines) and the six default
flips (commit after e31c029f). OPEN from the flips: (i) rt_split = true fatals when nx1 > RT_NNC = 72 on
the grey picket fence (all in-repo picket inputs are nx1 = 64) -> raise RT_NNC or make it a runtime size;
(ii) rt_split true vs false is NOT bitwise on the picket fence (1-KE 7 % on deep_hot_jupiter_rt_ideal_xe,
20 cycles, deterministic) contrary to the old comment -> a real discrepancy in the serial vs split grey
sweep, unexplained; only affects rt_ck = false runs. Both need a look before any grey dhj run.
