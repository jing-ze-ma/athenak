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
