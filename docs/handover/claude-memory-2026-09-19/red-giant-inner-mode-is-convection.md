---
name: red-giant-inner-mode-is-convection
description: RESOLVED 2026-09-13 -- the red-giant "inner radial mode" (e-fold ~3e5 s, i 5-40 above the wall, present with WB on/off/bounded) is PHYSICAL CONVECTION switching on at the base of the convection zone: the run's own Rosseland table gives nabla_rad/nabla_ad = 1.2 (i=5), 2.0 (i=7), 4.1 (i=10), 20 (i=20), 325 (i=40) so i>=5 MUST convect, and with mlt_alpha=0 nothing else carries L; the real-EOS superadiabaticity at the base built from 4e-5 (9.2e5) to 6.4e-4 (1.94e6), matching the measured rate; amplitudes 0.5-0.9 of (F/rho)^(1/3); i=5 already saturated by 2e6. The IC column (mlt_alpha_ic=3) is adiabatic to 1e-6 and motionless, so convection had to grow from the seed. The dynamic WB scheme only SHAPES it (C(v_r,T') 0.9 with it)
metadata:
  type: project
---

Trail: bench/RG_fofc_long2/mode (arms: only WB-off flattened it, over a window that only saw the WB-pattern decay),
mode2 (rmin/closure/plainsrc), bench/RG_fofc_long3 (wb_rmin=6e11: decay then REGROWTH at the same rate),
bench/RG_fofc_long3/analysis/{parta_realeos,partb}.txt (the real-EOS gradient with the eosgrad driver; the earlier
proxy verdict "8x too slow" used np.gradient on the 9e5 profile before the excess had built). Structure: peaks at
i 7-10 (r 1.6-1.8e11), wall cells pinned, smooth large-scale buoyant plumes. Remaining oddities: fluctuation
enthalpy flux 50-120x L/4pi r^2 at i 6-12 (T'/T ~2e-4, over-driven entropy perturbation?) and a slowly growing
(e-fold 8e5-1e6 s) angularly coherent odd-even MEAN flow in wall cells i 0-2 (ghost fill takes rho, e from the FIXED
initial column and mirrors only v_r) -- arms bench/RG_fofc_long3/mode: inner_open (Stein-Nordlund open bottom),
wallnoflux_off. **How to apply:** do not "fix" the mode; it is the convection zone turning on. Decide physics: run
with mlt_alpha>0 to carry the flux where the grid cannot, or accept resolved convection and its saturation. A
`wall_mirror` option (mirror the full wall state) is the missing numerical arm. wb_rmin/wb_rmax (hydro+MHD) and the
repaired static WB remain useful tools, not fixes for this. See [[red-giant-3e6-state]].

**Static WB tests (2026-09-13 ~10:30):** bench/RG_static (from scratch, wellbalance_static + static_reconst, NO wb_rmax
window) DIED at 5.63e5: dt collapse at a photosphere cell r 3.43e12 next to a block edge, T clamped at the table top
(2.629e15 code = 10^7.5 K), v_r -1.3e9; FOFC 1e6-5e6/interval from 3.5e5 -- the unbounded static background makes
the evolving photosphere/corona O(1) deviations. => static WB needs wb_rmax (+ wb_ramp) exactly like the dynamic
scheme; the base motivation is gone (the mode is physics). The restart arm (mode2/arms/static, also unbounded) was
left to finish for the record. RG_fofc_long3 (dynamic, wb_rmin=6e11) reached ~2.9e6 with KE_r 1.8e44, 3x long2's:
resolved base convection is MORE vigorous without the WB reconstruction there.
Static restart arm (mode2/arms/static, 1.0e6 -> 1.7e6, unbounded): the base mode GROWS x15.8 at i=20 (e-fold 1.76e5 s),
faster than the dynamic control -- consistent with physics that no WB variant suppresses; fofc ~0, tclamp ~100 (no corona
trouble on this restart, unlike the from-scratch run). Both 3e6 runs done: RG_fofc_long2 (dynamic WB) and RG_fofc_long3
(wb_rmin 6e11); comparison agent -> RG_fofc_long3/analysis.
Wall arms (RG_fofc_long3/mode, 1.4e6 -> 2.0e6): inner_open (no wall, Stein-Nordlund open bottom, flux 0) e-fold at
i=10 2.12e5 s vs control 2.19e5; wallnoflux_off 2.16e5 -> neither the wall nor its correction drives the base
convection. CLOSED as physics. All global RG jobs stopped by the user 2026-09-13 ~12:00; then RG_v4 launched
(bench/RG_v4, 11670310/1): from scratch, NO well-balanced scheme (user decision), binary 05130338 with the wall
energy fix, everything else as RG_fofc_long3; gates in NOTES.md.
