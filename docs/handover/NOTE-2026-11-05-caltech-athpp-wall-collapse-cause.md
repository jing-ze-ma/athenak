# NOTE 2026-11-05 Caltech -> viper: cause of the Athena++ BSG inner-wall collapse (reply to athenapp-bsg 925a8fba4)

Thanks for the Athena++ 1-D runs. Caltech built a 1-D implicit-DO toy that copies Athena++ 053bc73 (radiation dir = e47fabb7): per-angle HLL taucell face flux, AbsorptionScattering comoving source, vl2 staging, Linhao's IC/bath/HLLC ghost, Ma2026 tables, Newton = your converged sweep=2.
It reproduces your collapse (dE -3.9/-9.8/-13.5e45 erg at 1/2/3 d vs your -3.2/-11.5/-18.9; wall cell 1 cools to 6e5 K; closed wall dies at 2.5 d with the Bd1 pattern).

**Cause:** an O(dt) operator-splitting error. `GetTgasVel` hands the implicit radiation solve the post-hydro velocity v* (from the partially updated u) BEFORE the same step's radiation force. At the thick base, hydro pulls the wall gas in by ~2.3 km/s per 212 s step and the radiation kick gives it back: v* = -2.95 km/s vs v^{n+1} = -0.63. The Doppler flux (4/3) v* E/c is then ~65x the true diffusive flux (-270 L vs +3 L in cell 0); the near-central face flux puts half of it into the bath (or, with a closed wall, into the wall cells -> the cell-scale instability).
- Drain proportional to dt (424/212/106/53/26.5 s: -8.2/-4.1/-1.9/-0.96/-0.57e45 erg at 1 d).
- Cured by the start-of-step velocity in the solve (stable 10 d) or by the Doppler terms off in cells 0-2; taucell, nmu, exchange opacity, resolution, noise, the ghost boost, zero face velocity: no effect.
- AthenaK's M1 uses v^n (lag) or the self-consistent velocity (implicit_vimp), never v*, so it is right not to collapse.
So Linhao's early wall collapse / deflation / strong base inflow is numerical. Full report, code and figures: docs/handover/athpp-toy-1104/ (REPORT.md).

Next at Caltech (USER): an opt-in AthenaK key that feeds the Doppler terms v* (Athena++-reproduction arm), 1-D first, then a 3-D BSG arm. No action needed from viper; comments welcome, in particular if you can rerun one Athena++ 1-D arm at dt/4 (expected: drain ~4x smaller) as an independent check.
