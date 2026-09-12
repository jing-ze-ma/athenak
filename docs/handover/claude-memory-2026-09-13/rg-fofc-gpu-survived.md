---
name: rg-fofc-gpu-survived
description: RG_fofc on the viper GPU (2026-09-12, bench/RG_fofc, job 11625391) REACHED tlim 9e5 in 55 min with vceil=0 and FOFC on -- through the 5.7e5/5.87e5/6.13e5 death window with dt flat at 30.667 s; the shell launch DID happen (L_rad,out/L up to 117x at 5.9e5); one dt dip to 19.2 s at 7.0e5 (11/24 vertex columns under-dense at the JOIN, healed by 8e5); FOFC peak 62k flags/interval at 8.6e5 = 6.5e-5 of cells; state at 9e5 settling. FOFC replaced the ceiling on this realization
metadata:
  type: project
---

Config: cs 6x32x32, nx1 480, 24 MeshBlocks 480x16x16 (11% faster than 96 of 480x8x8), 1 node 2 MI300A,
109 ms/cycle, hydro/scratch_level=1 (LDS cap at nx1=480, 955c39df), binary build_hip_rg (da9734fd fixed
the sponge host-global capture). Input = orion's inputs/hydro/red_giant_fofc.athinput + those two lines.
Analysis in bench/RG_fofc/analysis (dt_table, counters_table, i8_out, shell.py/shell.txt, summary_extra).

Findings: dt min 19.2 s at 6.9989e5, otherwise 30.667 flat; zero "### dt COLLAPSE" blocks; eos_fail 0.
Floor hits (dfloor up to 1.9e6/interval) are in the CUBE-VERTEX CORNER GHOST cells: the quadratic
Lagrange corner extrapolation (10,-15,6) goes negative in the inflating layer above the join
(bvals_cc.cpp CsFillPanelCornersCC); active cells were clean in every dump except 16 cells at 7.0e5.
Hydro stencils never read that block ([[cs-wire-fill-wip]]); efloor_de total 1e12 erg = 1e-38 E_tot.
Mass above the join grew 7 decades 5e5->8.4e5 then declines; corona above 4e12 still loading at 9e5.
No classic vertex chimney (i 295-320 vertex columns are DENSER than bulk at 9e5).
max |v| (Cartesian norm) 3.3e6-1.26e7, always >= 4x under the retired 5e7 ceiling.
NOT determined: where the FOFC flags are (no flag field in dumps); a FOFC-off GPU control; anything
past 9e5 (rst every 1e5 to 8e5 in bench/RG_fofc/rst). The CPU twin bench/RG_fofc_cpu (11626142, apu host
cores) was launched to check the dfloor onset vs orion's "no event at 2e5" claim -- read its rg.log.
**Why:** the FOFC plan's purpose was to replace the velocity ceiling on the red giant; this is the first
survival of the death window without it. **How to apply:** compare with orion's own RG_fofc (CPU, same
input) when its links 2-3 have crossed 6.2e5; GPU and CPU are different chaotic realizations.

**2026-09-12 (user, from orion): orion's own RG_fofc (CPU, same input, links 196662-4) ALSO SURVIVED the
death window.** Two independent realizations (CPU orion, GPU viper) with vceil=0 + FOFC pass 5.7-6.13e5 ->
FOFC replaces the velocity ceiling on the red giant. Remaining: FOFC-off control on the same binary is
still absent (the deaths were older binaries with the ceiling); the ghost-corner dfloor onset vs orion.
**CPU twin READ 2026-09-12 (bench/RG_fofc_cpu, apu host cores, cancelled at c5543):** rg.log rows c4891/5217/5543
IDENTICAL to the GPU (dfloor 0/504/4898, efloor 76/2650/13698) and hst KE agrees to 4 digits at matched t
(2.4400e42 at 1.6e5, 2.3869e42 at 1.7e5): the red-giant pgen + FOFC run the same on MI300A as on CPU; the
ghost-corner floor onset is physics of the inflating join layer, not a GPU artefact. scratch_level 0 vs 1
bitwise on CPU (50 cycles). The GPU setup is validated for long runs.

**2026-09-12 ~08:00 CORNER LIMITER + LONG RUN.** Commit 2e3ae96e: sign-preserving fallback in FillPanelCornersCC/FC
(quadratic extrapolant replaced by the adjacent node only when the 3 nodes share a strict sign and q flips it).
Red giant restart from 8e5: dfloor 1.1e6 -> 320/interval, efloor 2.4e6 -> 3e5, history unchanged; raddiff cs
test unchanged to 4 digits (a node-range clamp broke it: Linf 0.036 -> 0.19, REJECTED). Remaining efloor hits are
the corona's own floor traffic (efloor_from_ekin), not the corner block. LONG RUN bench/RG_fofc_clamp SUBMITTED:
links 11633164/11633165 (apu, 24 blocks, 2 MI300A), restart from RG_fofc rst/rg.00009.rst (t=9e5), tlim 3.0e6
(~7 h at 109 ms/cycle). Gate: rg.log floors/fofc, dt, RG_fofc/analysis/shell.py + floormap.py (corona loading,
M(>join), vertex percentile ranks), the orion prod12 gate (open top + sponge off to 3e6).

**2026-09-12 ~15:00 LONG RUN DIED at t = 1.1192e6 (RG_fofc_clamp 11633164, killed 11633164/5).** dt 30.67 -> 1.22 in ONE
cycle (36598), then to 1e-26 by 1.124e6 spinning (NO time/dt_min in the red-giant input: ADD dt_min). The hydro
dt-cell diagnostic named (m,k,j,i)=(3,18,18,443) gid 15: r=3.906e12 (corona, above the 3.65e12 join), rho 1.2e-13,
T=1.3e17 K, v_r=9.8e8, cs 5.9e8; rad_cap_ang capping x_i up to 9e8 (cap 0.5) at i 422-443 of the same column;
rt_cell_report shows ei_u0 = -4.3e5 (negative conserved e_int, rt_use_cons clamp) at T=0.1 K cells nearby.
fofc rose 10 -> 3e6/interval between c32058 and c34667 (t ~1.06e6). = the corona thermal runaway at the loading front
(M(>4e12) was still climbing at 9e5), same family as the 1e16-K "vacuum cell" deaths in the red-giant notes.
Open: clamp or physics? -> RG_fofc_ctl (11639634, UNCLAMPED 955c39df from rst 9e5 to 1.3e6, dt_min 1e-3) and
RG_fofc_clamp/diag (11639635, clamped HEAD binary from rst 1.1e6, hydro_fofc map + hydro_w every 2e3 s to 1.125e6).
NOTE out.txt is binary-ish (rt_cell_report tags): grep needs -a.
**2026-09-12 ~16:00 restarts are NOT faithful (and the GPU IS deterministic).** determ job 11640095: two identical
restarts bitwise. But RG_fofc_clamp's own hst at t=1.11e6 vs the diag restart from rst 1.1e6: 1-mom 2.74903e36 vs
2.74902e36, KE 8.50297e41 vs 8.50296e41 after 326 cycles, and rg.log at cycle 36297: dfloor 6.8e6/efloor 3.7e6/
fofc 3.2e6 (continuous) vs 1.2e6/7.6e5/5.8e5 (restart) -- 5x fewer. Some evolving state is not in the rst
(candidates: the stateful open-top ghost fill, WB cache, RT guard state, FOFC-related arrays); the restart's
cleaner state delays the corona collapse from 1.119e6 to ~1.125e6 (diag hits dt 1.8 s at its 1.125e6 end).
Same pattern as the cs dhj deaths not reproducing from restarts (bench/cs_death_ab: all 3 arms survived).
Test running: RG_fofc_clamp/rstfaith (continuous vs restarted, bitwise cmp after the restart point).
**CORRECTION 2026-09-12 ~16:40: restarts ARE faithful.** rstfaith job 11640743 (same binary: continuous 400-cycle run
from rst 1.0e6 writing a mid rst; restart from it): hst BITWISE identical after the restart point. The earlier
"original vs diag restart differ in the 6th digit / 5x fewer floors" was ORIGINAL binary (build_hip_rg_clamp @2e3ae96e)
vs REBUILT binary (@7d87f3c5, + fofc_cnt kernel): hipcc kernel changes move results by ulps
([[hip-dualview-sync-idiom]]) and the corona runaway amplifies ulps within 300 cycles. Same for bench/cs_death_ab (HEAD
binary vs the production images). So: GPU deterministic, restarts exact, and the collapse TIME is chaotic at the ulp
level -- the collapses are real instabilities (red giant: corona thermal runaway at the loading front ~1.12e6; cs dhj:
unknown cell, rot 11-12 on two arms, none on cs_mhd_prod3 to rot 42). Audit sub-result: the only dynamics-affecting
state not in the rst is red_giant.cpp fmlt1d_ (MLT relaxed flux, re-seeded on restart, 1e4 s smoothing) -- did not
show in the bitwise test, so it is inert in this input. Open-top fills are memoryless.
**2026-09-12 ~18:00 CLAMP CLEARED, LONG RUN RELAUNCHED.** RG_fofc_ctl (unclamped) AND RG_fofc_clamp2 (clamped, rebuilt
binary) both survive 9e5 -> 1.3e6 with dt flat through the corona episode (FOFC peak 4.33e6/interval at cycle 35319 in
both, within 0.1%); only the first clamped build died there. => the episode is real, survival is ulp-chaotic, the
corner limiter is not implicated. RG_fofc_long (bench/RG_fofc_long, see NOTES.md for job ids) = clamped binary from
9e5, tlim 3e6, dt_min 1e-3 added to the input. Lesson: ALWAYS put time/dt_min in the red-giant input.
**2026-09-12 ~19:30 RG_fofc_long DIED at 1.7728e6** (dt_min FATAL, chain stopped; rst every 1e5 to 1.7e6 in bench/RG_fofc_long/rst).
Collapse cell (hydro dt diag): gid 15 (rank 1, m=3), (k,j,i)=(9,18,339), r=3.44e12 (photosphere), rho 5.9e-13, T=4.9e17 K,
cs 1.1e9, v_r=-8.4e6. The earlier death (old clamped build, 1.119e6) was gid 15 (18,18,443). With nghost 3 and 16-cell
blocks, j=18 is the LAST ACTIVE ROW = a PANEL EDGE and (18,18) the block corner = a CUBE VERTEX (2x2 blocks per panel).
=> seam/vertex T runaway near the photosphere, the [[red-giant-vertex-chimney]] family. Jobs: RG_fofc_ctl2 11649849
(unclamped 9e5 -> 1.85e6: does it die at the seam too?), deathmap 11649850 (deterministic restart from 1.7e6 with
hydro_w + hydro_fofc maps every 2e3 s incl. ghosts; analysis agent on it).
**CORRECTION 2 (2026-09-12 ~21:00), dump-level test bench/RG_fofc_long/bitwise (11650986):** two identical 300-cycle
restarts are BITWISE identical in the full hydro_w dump (GPU deterministic, confirmed properly). A restart from a
mid-run rst is NOT bitwise: after 150 cycles 1.2% of cells differ (36k of 2.95M), median rel 4e-6, max O(1), ALL in
the top ~30 radial rows (i 447-479, corona / open top); 98.8% of cells bitwise identical => not chaos (would spread),
not an ulp seed (median 1e-6): the running w0 differs from c2p(u0) that a restart recomputes (floors / RT clamp /
top-cell guard acting on w0 only?), or a top-boundary state. The earlier "faithful" claim rested on 6-digit hst over
200 cycles -- too weak. One-cycle test (rst1cyc 11651277: continuous vs restart dumped at the restart cycle and +1)
pending. The deathmap reproduction hit the same cycle/cell but dtold 17.7 vs 19.1 -> consistent with this.
**Restart bisection round 1 (bench/RG_fofc_long/bitwise/bisect, 11651679): NO single physics switch restores bitwise
continuation** -- base, rt_off, cond_off (rad_kappa_fac=0), implicit_off, fofc_off, open_wall, floors_off,
topclamp_off, usecons_off, wbcache_off all DIFFER after one cycle. One-cycle facts: at the restart cycle ACTIVE cells
are bitwise identical; +1 cycle: eint differs in i 384-479, dens in i 444-479, velx everywhere (median 3e-5).
Round 2 (bisect2): ghost-inclusive comparison at the restart cycle, WB off, etotgrav off, reconstruct=dc, 1 rank,
plus a code read of the first-cycle-after-restart path (Driver::Initialize BC/prolongation/c2p ordering vs the
stage task list; first dt from the file?). Suspect: the ghost/seam state or an init-path asymmetry, not physics.
**Restart bisection round 2 (bisect2, 11652242) FOUND IT: at the restart cycle ACTIVE cells are bitwise identical but
the OUTER RADIAL GHOST layer (i 483-485, the open top) differs in ALL 18432 cells by rel 5 -> 1e6 -> inf, angular face
ghosts near the top by 1e-3; one cycle later the top active rows differ. Same in every arm (wb_off, etotgrav_off,
recon_dc, onerank, wall BC smaller). => the restart's first physical-BC fill (Driver::Initialize chain
ApplyPhysicalBCs -> ConToPrim, same order as a stage) computes the open-top ghosts from a different/invalid
input (w0 stale-by-one-stage in a run vs never computed at restart init?) than the running run; the rst file DOES
carry the running ghosts and the init chain overwrites them. Agent fixing in red_giant.cpp (BC to read u0 or skip
the overwrite). Same lag likely in deep_hot_jupiter_rt.cpp's user BCs -> every chained dhj run restarts non-bitwise.
