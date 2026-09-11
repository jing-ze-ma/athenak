---
name: red-giant-stale-primitive-audit
description: "STALE-PRIMITIVE AUDIT 09-11 02:40 (read-only, all of red_giant.cpp user srcs/BCs, two_stream_rt, srcterms, diffusion, eos): the RT bug class (split operator bounds/relaxes against w0 while writing u0) recurs at 6 more sites; LIVE in the current config: rt_use_cons default false (arm it), RedGiantBC outer/inner ghost fill seeded from the stale w0 of the active cell (B, at exactly the faces where runs die), rg_wallflux eos.Pressure on unguarded ei (C). LATENT (not enabled in the runs): rg_relax (grey relax, same class), rg_co5A/B/D open-inner BC (stale targets + unguarded ei + ORTHOGONAL KE on the cs), MLT 0.1*e/dt caps at 2305/2369, hydro_tasks perturbation-reconstruction ordering (conduction/viscosity/FOFC read the perturbation w0)"
metadata:
  type: project
---
Principle: explicit RK sources read w0 (stage input) and add beta_dt*S(w0) to u0 = CORRECT
(gravity/WB 1918, MLT flux, explicit conduction fluxes, radcapx 678, srcterms, turb driver);
split/implicit operators (relax to equilibrium, guards, energy-fraction limiters, floors,
damping) must bound the CURRENT u0 (EintFromCons) = ImplicitRadialUpdate is the reference.
Ordering correction: user_srcs is called INSIDE HydroSrcTerms (hydro_tasks.cpp:517), i.e.
BEFORE ImplicitConduction, and within RedGiantGravity the order is FaceBudget -> gravity/WB
-> MLT -> open-inner BC -> wall no-flux -> sponge -> RT -> grey_relax; later w0 readers are
stale by all earlier u0 edits too. rg_mlt_shell uses the PREVIOUS call's rt_face_flux (lag).
Ranked list (file:line, red_giant.cpp unless noted):
1 two_stream_rt.hpp:332 rt_use_cons=false default -> set true in inputs / flip default.
2 RedGiantBC 2953-3051: ghosts built from w0 of the active cell after split ops moved u0
  -> derive d, e, v from u0 (EintFromCons + CsKinetic), keep the open-ghost guard. LIVE.
3 rg_wallflux 2692-2696: eos.Pressure(d, ei) with no ei>0 guard. LIVE (inner_bc=wall).
4 rg_co50 2451, rg_co5B 2524-2527 (also orthogonal KE 0.5*r*(v1^2+v2^2+v3^2) instead of
  CsKinetic -> wrong cross term on the cs), rg_co5D 2600 (same KE), rg_co5A 2491 (meanP,
  srho0 targets from stale w0). Latent (inner_bc=wall in the runs).
5 rg_relax 2823-2839: exponential relaxation with e from w0. Latent (rt_grey on).
6 MLT caps 2305/2369 fmax = 0.1*e_w0/dt. Latent (mlt_alpha=0).
7 hydro_tasks.cpp:296 vs 485 (use_wellbalance_static_reconst_perturb): w0 holds the
  perturbation during conduction/viscosity/FOFC. Latent. Mirror in mhd_tasks 191/290.
Clean: atm_column.hpp, correlated_k.hpp, viscosity, resistivity, hydro_update, general_hyd
floors (u-based; wtemp guess only). See [[red-giant-runaway-source-rt-stale-w0]].
