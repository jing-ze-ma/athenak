---
name: sparc-sponge-campaign-0925
description: 09-25 plan -- SPARC-like C32 hydro fresh start 300 rot (nx1 ~60-72, 3/H below 1e-6 bar, 8-coef) + 3 sponge-test arms 100 rot; dir /viper/ptmp2/jinma/sparc_0925
metadata:
  type: project
---
User 09-25: (1) "SPARC-like" run: C32 hydro, fresh start from the analytic IC, radial grid 3 cells per median H for
p > 1e-6 bar (free above), 8-coefficient stretch, T4 + c2 + ck_impl_every = 4, tlim 300 rot -> CHANGED 09-25: 1 GPU and 100 rot like the tests; the arm chosen as final baseline is extended from its restarts; shorter links
if the queue waits. (2) Sponge tests at the same grid, 100 rot each: T1 bottom sponge off, T2 top sponge off, T3 both
off (switches problem/sponge_top / sponge_bottom, branch dhj-sponge-switch, default on = bitwise). UM-like uniform-r
arm (nx1 168) designed but NOT run. Radial A/B a128 done, a320 queued (start est. 09-25 08:13) -> decides the
production spin-up grid (256 8-coef recommended, see [[mhd-relax-hydro-spinup-ok]]).
Sponges today: top = all 3 momenta, 0 at 1e-6 bar -> 1/1000 s at <= 1e-7 bar; bottom = horizontal momenta only,
0 at 50 bar -> 1/1000 s at >= 100 bar (very strong vs literature basal drag 1e5-1e7 s); plus the 5-rotation spin-up
Rayleigh drag. BCs: top isothermal hydrostatic ghosts + zero-gradient velocity; bottom ghosts held at the IC state.
Binary for the whole campaign fixed (commit + md5 in the run dir). Smoke on apudev before each chain.
**Update 09-25 (user):** the four runs WAIT for a new initial condition: 1-D ck radiative-convective-equilibrium profile,
GLOBAL-MEAN irradiation (f = 1/4), exact EOS adiabat below the RCB (replaces the substellar picket-fence + 0.9 nabla_ad IC),
via new problem/ic_profile (branch dhj-ic-profile, /viper/ptmp2/jinma/ckrce_0925). Deep-convection diagnostic in
/viper/ptmp2/jinma/deepconv_0925. Sponge switches merged 6375568c.
**HOLD (user 09-25, later):** initial-condition (ck RCE) work and the sponge runs are ON HOLD; ck_nquad=2 validation continues. IC agent stopped at a clean point, state in /viper/ptmp2/jinma/ckrce_0925/README.md.
IC resume traps (ckrce_0925/README.md): ic_profile table must reach >= 250 bar (get_init_Tp_host assumes uniform log p and reads past the end); in blend mode the ck cut = first face with w < 1, pcut ignored.

**User 09-25: the four sponge tests run on WASP-121b** (wait for /viper/ptmp2/jinma/wasp121_0925 inputs + the rebuilt binary; then smoke, then submit).
- a320 resubmitted 09-25 08:10 as 1 node (2 GPUs) x 4 h, WALL 03:40, link 1 11971805 + link 2 11971806 (afterany, restarts from newest rst); same binary as a128 on purpose (md5 5ba3c5df).

**User 09-25: only p > 1e-6 bar data matter; the model top may be lowered if the flow there is unaffected.** Test: two extra arms (base sponges) with p_top 1e-8 and 1e-7 bar vs baseline 1e-9, compared on the flow at p > 1e-6 bar. Note the top sponge window is hard-coded 1e-7..1e-6 bar.
Top-boundary arm SKIPPED (user 09-25): lowering the top to 1e-8/1e-7 bar saves only 2 cells (nx1 74->72) and no dt; keep p_top 1e-9 unless the top causes floor/convergence trouble. Inputs sparc_w121_top8/7 kept in wasp121_0925.
- User 09-25 ~22:30: sponge runs use the TESTED binary (sparc_0925/athena.gpu md5 e3a8442e, without ck-fast2 levers); ck_dif_dtau=10 + ck_beam_par go into later production runs.

**SUBMITTED 09-25 ~22:45: WASP-121b 1x sponge runs** (PLANET=w121x1, binary e3a8442e, 1 GPU per arm, TROT 100): base+nobot 11980312 -> 11980313 (afterany), notop+none 11980314 -> 11980315; arm dirs /viper/ptmp2/jinma/sparc_w121x1/<arm>, logs there. Smoke 11980271: FATAL 0; ck 9/75 nonconv with top sponge (0 in restart), 0 without; notop/none dt 15.5 -> 4.4 s in 300 cycles; floor events ~millions (night top). Check after a few rotations: T-floor cells below 1e-6 bar per arm.
- 23:17: sponge runs STARTED early on apudev: chained 14-min chunks (WALL 12 min), base+nobot 11980812->13->14, notop+none 11980815->16->17; apu links 11980312/11980314 now afterany the last chunk and continue from the newest rst (06:48 est.).
- 23:25 user: don't use apudev for the sponge runs; cancelled the pending apudev chunks (11980813/14/16/17); the two running first chunks finish; apu links now afterany 11980812 / 11980815.
- 23:28 user: 3 apudev chunks per pair is fine, not more (re-added 2 per pair after the running ones).

09-26: arm none (no sponges) dt-collapsed at cycle 3777 (t 1.99e4 s = 0.18 rot, hydro dt 0.5 -> 0.03 s); user stopped it (none/STOP). notop dt 8.4 s vs base 13.6 at 1.6 rot. base/nobot healthy (5.4 / 3.1 rot).

**T-floor check 09-26** (/viper/ptmp2/jinma/sponge_chk_0926, tfloor_chk.py on the restarts, tfloor 200 K): at p > 1e-6 bar every arm has 1-11 floor cells (<= 3e-5 of cells) at 1-20 microbar, night side, |lat| 15-26 deg; the count is about the same in base/nobot/notop. The top sponge only buys dt (13.6 vs 7.5 s; without it |v| is 19-51 km/s at 1e-11 bar night). The bottom sponge changes nothing so far (<= 8 rot). Floor events grow: 2e7 -> 1.3e8 per rot (notop 2.4e8 at rot 4). The floor cells above 1e-6 bar (0.5 %) never set dt.

**DECISION (user 09-26): keep both sponges. Recheck the bottom sponge (base vs nobot, deep flow and T at 1-100 bar) at rot 20;** rerun sponge_chk_0926/tfloor_chk.py as well.

**09-26 ~02:30: nobot and notop CRASHED (FATAL dt collapse, STOP written).** nobot at rot 6.03 (t 6.649e5, cell gid 23 i 50-54); notop at rot 5.59 (t 6.160e5, gid 0 i 56). In both, the cell is at the density floor 1e-16 and the pressure floor 1e-5 (cgs), T 1e11-1e13 K, |v| 1e7-3e11 cm/s: an evacuated upper-atmosphere cell. Only base survives (rot 18.6, dt 13.6). So the rot-20 base-vs-nobot recheck is impossible as planned.
- PLANNED (user 09-26): when the crash diagnosis (/viper/ptmp2/jinma/sponge_crash_0926/RESUME.md) reports, retry nobot from its last restart before the crash, with the fix it suggests if any. Pair it on the idle second GPU of the base chain (the next base link 11980313), or give it its own 2-arm job.
**Crash diagnosis 09-26** (/viper/ptmp2/jinma/sponge_crash_0926/RESUME.md):
- Both crashes are on cubed-sphere face-edge SEAMS, on the night side, locally at 1e-7..4e-7 bar (inside the top-sponge window). 6 of 7 dt-collapse sites are on face edges (12 % of cells).
- The reruns reproduce each crash exactly.
- nobot: one seam row pinned at the 200 K T-floor next to 1000-7000 K neighbours; the along-seam flow reaches M ~55, then the cell blows up.
- notop: a supersonic night downdraft (M 10-20) drains the seam column, with spurious heating (E - KE round-off).
- The end of both: the legacy density floor raises rho but keeps momentum and energy, giving v 3e11 cm/s, T 1e11+ K, dt collapse.
- base has the same precursor state and no collapse yet. The bottom sponge is NOT shown causal (one event).
- Suggested fix, being gated: <hydro>/dfloor_keep_velocity, dfloor_keep_temperature, vceil ~1e7.
- Open: why a seam row stays at the T-floor (ck seam face mixing or the seam ghost exchange?).
- Floor-fix gate 09-26 (job 11982198): with <hydro>/dfloor_keep_velocity, dfloor_keep_temperature and vceil=1e7, both crashes are gone (nobot and notop ran past them, no warnings). vceil fired 355 / 1946 times near the old crash; floor counts unchanged. Stability only; accuracy not shown. User 09-26: LEAVE the decisions (continue nobot, production floor keys, base) for after the handover.
- **09-26 ~03:00: base (both sponges) ALSO CRASHED** at rot 27.04 (FATAL dt collapse cycle 221708, t 2.978282e6, gid 1, cells (1,10,17,60) / (1,2,17,49); 7 collapse warnings before it). ALL arms of the campaign have now died. Sponges only delay it; the seam/floor defect is general. Seam agent informed.
