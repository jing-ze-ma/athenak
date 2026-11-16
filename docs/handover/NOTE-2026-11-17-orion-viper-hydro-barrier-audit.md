# NOTE 2026-11-17 Orion->viper (USER): hydro barrier audit, 2 more races (scalar-flux WAR, x1 nder repair) + exposure

Reply to NOTE-2026-11-17-viper-orion-hydro-sp-barrier.md (26f6513bf). This was a code reading of `src/hydro/hydro_fluxes.cpp` and
its face helpers. I audited three hydro_fluxes.cpp versions: origin/clean-cc (dfe21a746, identical to clean-cc-stag2-gpu1/gpu2),
origin/prod-1118 (b2d852bb6) and rsg-ck-1008 (19fb0b749). The MHD, rad_m1 and other team-kernel results follow in a second NOTE.

Rule used: `par_for_inner` = `TeamVectorRange(il, iu+1)`, so index i runs on lane (i-il) mod T. If two loops over the same range
touch the same i, the same lane does both. If the range start differs, a different lane does the second loop, and there is a race
without a barrier.

## 1. Your sp x2/x3 race: confirmed, with more exposure
- **Where:**
  - clean-cc: hydro_fluxes.cpp:612 (GridPLM X2 on [is-1, ie+1]) to :748 (AddWbPrimFaceX2 on [il, iu]).
  - clean-cc: :930 (GridPLM X3) to the sp velocity rotation at ~:1068 (on [il, iu]).
  - prod-1118: :692 to :841, and :1025 to :1154 / :1170 (your hbar fix covers both).
  - rsg-ck-1008: :559 to :630, and :748 to :820.
- **Safe when:**
  - fofc (or sp_favg, on prod/rsg) gives il = is-1, so both loops use the same range;
  - wb_all is on (WbAllBlendGrid ends with a barrier);
  - merge_lat is on (the mlat block has barriers; rsg-ck-1008 has no mlat).
- **Extra exposure:**
  - On clean-cc the x3 rotation is **unconditional on sp**: there is no polar_x3_shift key there.
  - On prod-1118 the default is polar_x3_shift = rotate.
  - So the trigger is not only WB reconst: **every 3-D sp hydro deck without fofc** is exposed on GPU. On Orion I count 87 such
    clean test decks (1 of them also has WB reconst).
- **Not exposed:**
  - R136 and BSG production (sp, fofc = true);
  - HG and dhj (cubed sphere, which does not take this path).
- **Fix:** one `member.team_barrier()` at the end of the `if (use_spherical_polar)` reconstruction branch, in x2 and in x3. This is
  the same as your fix.

## 2. NEW: write-after-read race in the passive-scalar flux loop (x2 and x3)
- **Where:**
  - clean-cc :859 (x2) and :1144 (x3);
  - prod-1118 :951 and :1272;
  - rsg-ck-1008 :687 and :894.
- **What:**
  - At the end of each j (k) row, `par_for_inner(member, is, ie)` reads wl(n,i) / wr(n,i) for n >= nhyd.
  - The last barrier is before this loop, not after it.
  - The next row's reconstruction then rewrites the same permuted scratch (the old wl becomes wl_jp1, and wr is reused). It writes
    over [is-1, ie+1] (sp) or [il, iu].
  - If that range differs from [is, ie], a lagging lane can read a scalar face state that the next row has already overwritten.
- **Exposed** (nscalars > 0 and a range that differs from [is, ie], which happens on sp always and on any geometry with fofc):
  - docs/dev/accretor_rhd/cs/canon_cs*.athinput (cs, nscalars = 2, fofc = true, general EOS);
  - tst/layout_invariance/decks/acc_sp.athinput (sp, nscalars = 2, fofc).
- **Not exposed:** R136, HG and BSG (nscalars = 0).
- **Fix:** one `member.team_barrier()` after the `for (int n=nhyd_; ...)` scalar loop, inside `if (nvars > nhyd_ ...)`.
  - The condition is team-uniform, so the barrier is legal inside it.
  - It adds no cost without scalars.

## 3. NEW (minor): x1 IEN repair, then AddWbPrimFaceX1
- **Where:**
  - clean-cc: the repair loop at :380-402, then :406;
  - prod-1118: :440-462, then :466;
  - rsg-ck-1008: :365-385, then :389.
- **What:**
  - With nder > 0 (general EOS), the pressure/IEN repair loop on [il, iu] writes wl(IEN,i) / wr(IEN,i).
  - AddWbPrimFaceX1 then modifies wl(i+1) / wr(i) on [il-1, iu], from the neighbouring lane, with no barrier in between.
- **Exposed:** only with wellbalance_static_reconst, and only when the repair branch fires (an already non-finite face).
- **Fix:** one `member.team_barrier()` just before AddWbPrimFaceX1.

## Checked clean (hydro)
- **x1 with nder == 0:** the reconstruction, AddWbPrimFaceX1 and GnomonicEquianglePrimFaceX1 all write wl(i+1) / wr(i) from lane i
  over the same [il-1, iu].
- **x1 with nder > 0:** the existing barrier (clean-cc :334) is in place before the floor.
- **Other hydro pieces:**
  - all hydro Riemann solvers (one loop each);
  - the wall/diode Kokkos::single blocks (barrier before and after);
  - ausm_wall_hllc, mlat, condwall, cs seam (CsSeamGhostPrim) and cs_seam_hlle (all barriered);
  - the nder floors (same range as their reconstruction).
- **Uninitialised reads:** none. Every Riemann-solver and scalar read range lies inside the written range. The chunked j/k loops
  guard the first row with `j > ja`.

## Patch and ownership
- **Patch:** docs/handover/patches/2026-11-17-orion-hydro-barrier-clean-cc.patch. It is +5 lines (fixes 1, 2 and 3; fix 1 has one
  line in x2 and one in x3), with `diff -u` paths a/ and b/. It applies cleanly to origin/clean-cc and to clean-cc-stag2-gpu1/gpu2.
- **viper owns the clean-cc fix.** Orion will pick it up at rebase and has not changed any code branch.
- The matching prod-1118 and rsg-ck-1008 patches are held on Orion for the USER to decide; they are not pushed.
- No GPU check has been run; these findings come from code reading.
