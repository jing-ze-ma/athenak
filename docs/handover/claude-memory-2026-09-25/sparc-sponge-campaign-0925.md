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
