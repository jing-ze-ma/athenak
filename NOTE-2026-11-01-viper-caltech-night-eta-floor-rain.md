# NOTE 2026-11-01 viper -> Caltech (USER): nightside eta is probably 80-800x too high (Saha x_e); s10 floor rain confirmed; ambipolar + x_e tracer coming

Reply to NOTE-2026-11-01-caltech-viper-mhd-floor-rain (8589f89d4), plus a new result that bears on r1m and on any implicit-resistivity work.

## 1. Nightside x_e is far above Saha, so the night eta_O is far too high
Offline estimate on our s10 (10 G, max_eta 1e15) rot-318 dump: backward streamlines from the night side to the day side, then forward ion chemistry along each path (Na, K, Ca, Al, Mg, Fe, Si, H; Verner radiative recombination, three-body, H- sink, condensation). Viper: /viper/ptmp2/jinma/nightxe_1101/RESULTS.md.

Median F = x_e(advected) / x_e(code Saha), zenith angle > 120 deg:

| log p [bar] | -7..-6 | -6..-5 | -5..-4 | -4..-3 | -3..-2 | -2..-1 |
|---|---|---|---|---|---|---|
| F, base (on isobars) | 3 | 12 | 315 | 780 | 83 | 1.2 |
| F, recombination x10 | 2 | 3 | 45 | 87 | 14 | 1.0 |
| F, with vertical motion | 106 | 248 | 523 | 896 | 141 | 1.3 |

- **The carriers are Na+ and K+** (60-90 % of night electrons). The other metals condense. The ions keep their dayside charge for the 8-21 h night transit, because radiative recombination of closed-shell ions is slow. Three-body recombination is negligible.
- The morning terminator re-ionizes at once (F 1.0-1.4). The dayside is in Saha (Damköhler 3-2e4).
- **Consequences:**
  - Night eta_O at 1e-5..1e-2 bar drops by F.
  - Night Rm goes from 0.06-6.5 to 9-2400 (median, length H_p).
  - The binding resistive dt limit in r1m (night, 1e-3..3e-2 bar) is therefore set by an eta that is probably 10-800x too high. Most of the STS stage count is spent on that layer.
  - **Please weigh this before investing further in implicit Ohmic resistivity.** With a realistic x_e the stiff layer mostly goes away.
- **Uncertainties:**
  - Missing channels (molecular-ion clustering + dissociative recombination, grains) could raise the effective recombination rate by an unknown factor.
  - F scales roughly as 1/alpha, so x10 still leaves F 14-87 at 1e-5..1e-2 bar.
  - Steady flow from one dump.
- **Hall:** with Saha x_e the night Hall/advection ratio is 25-220 (median). With the advected x_e it is 0.05-0.27, but 13-30 % of cells stay above 1, and above 1 again with alpha x10. Not implemented; marginal.

## 2. In progress on viper (the user's go, 11-01)
- **Advected x_e tracer.** One passive scalar, z = x_e, relaxing to the local Saha value with d n_e/dt = alpha(T)(n_e,eq^2 - n_e^2) (implicit, pointwise). It replaces x_e in eta_O and in the ambipolar rho_i.
  - Planned key `<mhd>/xe_tracer`, default off = bitwise; it will start from restarts without the scalar.
  - We will validate on a W121 MHD restart against the offline F and time it on MI300A. We will send a NOTE before it goes to fin51-r-fix.
- **Ambipolar diffusion** (`<mhd>/ambipolar`, default off = bitwise). Done in Cartesian, sp and cs. It is being merged into fin51-r-fix now, on top of your rkg_lean / rkg_dt_safety (899153cb), with key-off bitwise gates on GPU.
  - Our offline estimate: negligible on the dayside and terminator. AD/advection ~1 at night at 1e-6..2e-4 bar, but that uses Saha x_e; the tracer will reduce it.
  - On a W121 restart, ambipolar = eos lowered the explicit-resistivity mean dt from 0.82 to 0.63 s; RKG barely changed.

## 3. Floor rain in our s10: confirmed
From s10_e15 dhj.user.hst (Mfloor column, rot 300-320):
- Mfloor = 3.4e15 g/s (median, rot 300-305), rising to 4.1e15 g/s (rot 310-320).
- That is 4.5e20 g per rotation, ~16 % of the mass above 1 mbar (4 pi R^2 p/g ~ 2.8e21 g), and more than the mass above 0.1 mbar.
- It matches your estimate (2e15 g/s at r1m-smoke B). It was not flagged anywhere on our side.
- **Implication for the WASP-121b Fe-wind number** (NOTE ...-w121-fe-wind-vs-seidel, 5b149b42a). In s10, Fe forms at ~0.1 mbar, inside the layer this rain reaches. Its Kp-Vsys match (-4.26 vs -4.75 +- 0.30 km/s) is now caveat 7: infalling floor mass at the Fe level, effect unknown.
- We support valfven_noinfall_fac. Please send the fac 4/10 results. The stage-1 B calibration (3/10/30 G) should run with it.

## 4. Nothing needed from you now
The nightside-eta point is for planning (r1m dt, implicit resistivity). The tracer NOTE follows.
