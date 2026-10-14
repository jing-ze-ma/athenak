# NOTE Caltech -> Raven: He giant FIN4 domain too small for the tlim-1e9 extension. Please plan a larger-domain run (USER 10-14)

USER (10-14): "does the he giant run have enough r max", then on the options "yes do both, ask raven".

## Where it stands (Caltech wide4/run, page https://claude.ai/artifact/9g8xcCyBCjdzwkBbALaiqt)
- Domain: x1min 3 Rsun, x1max 200 Rsun (outflow top), nx1 897, sponge r0 180 Rsun (he_sponge_r0 1.25226e13, rate 1e-3, dmax 1e-17).
- R_ph (shell-mean tau_R = 2/3): 61.9 Rsun at t=0, then 74.6 / 82.4 / 88.2 / 88.5 at 10.75 / 12.75 / 15.75 / 16.25 d. The outer surface r(<rho> = 1e-12) is 98 Rsun at 16.25 d. vmax is 39 km/s at 107 Rsun. Expansion is about 1.2 Rsun/d and slowing.
- Extrapolation: the surface reaches the sponge (180 Rsun) at roughly 85-95 d. tlim is now 1e9 with 7-d links (about 240 d of star time per link: 4513229, 4513230 after 4489939, which stops at 30 d), so the first extension link would cross it.

## Done on Caltech (option 1)
- **Guard job:** hgGuard 4514966 (CPU, expansion, 7 d; follow-on 4514967 afterany) checks hourly. If r(<rho> = 1e-12) >= 150 Rsun or R_ph >= 150 Rsun, it writes wide4/run/STOP and stops every hgF4w job; the chain's link script honours STOP. Script: hegiant_1014/rmax_guard.sh (RLIM env, default 150).
- First line: t = 16.25 d, r12 = 98.1, R_ph = 88.5 Rsun.

## Request (option 2): plan a larger-domain He giant run
Please write a TASK for Caltech (or decide where it runs) for a FIN4 He giant with room for the long-term expansion. Points to decide:
1. **x1max:** e.g. 400-600 Rsun, or what a projection of R_ph(t) needs for about 1 yr of star time.
2. **Grid:** nx1 and the radial ratio, keeping the current resolution inside about 100 Rsun (the photosphere and convection zone). Also the cost on 4 H200: 897 cells now, about 0.5 s/cycle, dt about 19 s.
3. **Sponge:** r0 = 0.9 x1max per the E1 rule; check that the outflow BC and floors are consistent with the new top.
4. **Start:** fresh from the MESA IC, or remapped from a Caltech rst (e.g. near the guard stop) onto the larger grid. A restart cannot change the mesh, so a remap tool would be needed. Which do you prefer?
5. Gate: the same CPU/GPU checks as FIN4, and tlim 1e9 with 7-d links per the long-tlim rule.
Caltech keeps the current chain running until the guard stops it.
