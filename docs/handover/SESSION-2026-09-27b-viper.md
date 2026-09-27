# Viper session state 2026-09-27 evening (supersedes SESSION-2026-09-27-viper.md) — START HERE on viper

Read this, then `MEMORY.md`, then `HANDOVER-2026-09-26.md` for the full project background.
Always `git fetch fork` before any push (the Caltech session pushes to the same `rt-integration`). Never force.

## 0. First actions, in order

1-3. DONE 09-27 evening: the combined merge (4 branches + wall_closed default true + fac33815 MHD null-phydro fix
   + 2bee21a2 thin_relax ctr_mem saved in the rst + merge of Caltech 0da9a439) passed every combo gate and is
   PUSHED: fork/rt-integration = c32d8bf1. Gate table in /viper/ptmp2/jinma/combo_0927/RESUME.md (MHD growth
   rate cfl 0.3 = 0.9: M0 1.77e-4, MW 1.43e-3; dhj closed wall Etot_bot = Lrad_bot every row; tst all pass).
4. **Restart WASP-121b 1x and 10x** (user wants it; confirm launch with the user): from
   `/viper/ptmp2/jinma/w121prod_0927/w1x/rst/dhj.00139.rst` and `w10x/rst/dhj.00131.rst`, with a binary built
   from the pushed rt-integration (closed wall now default; ck_nquad 2 / T4 defaults change nothing there since
   the inputs name them — verify with `-n` parameter dump that the effective keys equal the old production
   keys except wall_closed). Delete the STOP files in w1x/ and w10x/ only at launch. `run.sub` +
   SUBMITTED.txt there. Goal still rot 300 before any remap. 2 GPUs per node, HSA env vars.
   Watch: 100-bar T drift (the wall was NOT its main cause: +5.4 K vs +4.5 K per 3 rot), 1-10 bar cooling
   (20-30 K / 3 rot in both arms), Etot_bot = Lrad_bot, Lir_top / (Lsw_abs + Lrad_bot) -> 1.
5. Cleanup (ask the user): `/viper/ptmp2/jinma/wallE_0927/build_gpu`, `build_dbg`, the failed first run dirs
   (A-G, T1-T9); `/viper/ptmp2/jinma/defaults_0927/{src_base,src_new,bin}`; old `kedt_0926/src_k*` builds.

## 1. Done this session (all verified by the main session against the run output)

| item | result | where |
|---|---|---|
| dhj flux history | merged 392f1e50; 8 cols in dhj.user.hst; Lrad_bot = sigma T_int^4 4 pi r_in^2 (1.002) | fluxhst_0927 |
| ke-dt coupling fix | merged + pushed 86d7c115; wedge gamma cfl-independent (-2.10e-5 vs -2.12e-5) | kedt_merge_0927 |
| ke-dt keys default (wb_arad only) | branch b5c43608; wb_arad simple tests: order 1.2 -> 2.1, v err at cfl 0.9 10-13x smaller | kedt_merge_0927/tb |
| inner-wall energy leak | ix1_bc=user never writes the fluid ghosts -> frozen reservoir; net 1.28e30 erg/s (137 L_int) at zero net mass flux; wall_closed -> Etot_bot = Lrad_bot, Mdot 0, box dE/dt 1.71e30 -> 3.9e29 | wallE_0927 |
| starlight columns | Lsw_in 6.33e30: refl 16.5 %, abs 43.2 %, 40.3 % transmitted through the thin top shell (physical); ck r^2 NOT missing | wallE_0927 |
| recommended switches default-on | branch 458160fb; full list + UNSURE in `defaults_0927/INVENTORY.md`; `docs/dev/default_flips_0927.md` | defaults_0927 |
| M1/VET on MHD | branch f1d3491a; gates (i)-(vi) pass; MHD growth rate cfl 0.3 = 0.9 (1.91e-4 B=0, 3.81e-4 77 G) | m1mhd_0927 |
| srclim tst fix | pushed 71332051 (reads hydro.hst; flux_hst adds user.hst) | — |

## 2. Open decisions / threads for the user

- **defaults UNSURE list** (not flipped): ck_implicit itself; c2 levers (ck_impl_once, jreuse 0.2, pred) and
  ck_impl_every 4 (not round-off neutral); ck_impl_floorbound, kkt_demax; ck_pcut_bar; grav_point_mass,
  rot_potential; f_source=wb; implicit_solver=bicgstab, vet_mb_tblock (need tests); ck_nquad 2 on red_giant /
  box_convection; marshak_face / vet_col_order2 on Cartesian; sponge_bottom.
- **m1-mhd remaining**: cs + MHD + M1 (FATAL today), box rt_two_stream / wb_arad_force under MHD, other rad_m1
  tests. Pre-existing bug (hydro too): box_convection wall_noflux subtracts rho*phicc_wb (effective potential)
  from an etotgrav energy holding the true potential.
- **dhj deep drift** after the wall fix: still unexplained (see step 3).
- Remap nx1 76 -> 256 after rot 300 (radab_0924); 10x ck denser Jacobian (RESUME-ck-newton10x-b.md).

## 3. Rules learned this session (in memory)

- Simple tests must exercise the new code path (keys that are no-ops do not count).
- Gate the COMBINATION of merged branches before pushing.
- Never `pkill -f` a pattern that matches your own shell; never amend in a checkout a worker also commits to.
