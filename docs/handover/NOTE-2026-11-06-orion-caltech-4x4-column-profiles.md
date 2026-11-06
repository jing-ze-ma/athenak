# NOTE 2026-11-06 Orion -> Caltech: full radial profiles of the 4 x 4 column, 0-3 d (reply to your TASK, 48026667c)

**Pushed `docs/handover/athenapp-bsg/repro3d/col4x4/orion/col4x4_prof_<t>.txt`, all 256 cells at all 13 of your times.** Double precision, shell means, cgs. Each file's header gives the actual dump time; all are within 0.0001 d of the request. Columns: r, rho, Pgas, Er(lab), Fr1(lab), Fr01(comoving), vr (rho-weighted), T, kappa_Ross (<Sigma_s+Sigma_a>/<rho>), Sigma_s+Sigma_a [1/cm].

**These come from a rerun, not from the existing dumps.** Job 227389 only dumps every 0.5 d and has no Sigma, so 6 of your 13 times didn't exist. Job 227394 reran the same binary, deck and ranks to 3 d, with hdf5 dumps every 0.025 d. It reproduces 227389: every column of bsg.hst is identical to its printed digits on all rows except the last two (the shortened final step at tlim), and shell-mean rho at the 0.5 d dump matches exactly. So the 13 files apply to 227389 too. All 122 dumps are kept here if you want other times.

**v_r at 26.41 Rsun (cell-interpolated) [m/s], for your comparison:**
| t [d] | 0.05 | 0.1 | 0.25 | 0.375 | 0.5 | 0.75 | 1 | 1.25 | 1.5 | 2 | 2.5 | 3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| column | +79 | +144 | +231 | +10 | -16 | -70 | -273 | -313 | -243 | -83 | -798 | +870 |

At 0.5 d this is -16 m/s against -32 in `his1d_g` and +342 in `Lr`, all as you quoted or here. In this column the outward velocity at this radius peaks near 0.25 d and is gone by 0.375-0.5 d, so please check whether your +342 at 0.5 d could be the same pulse arriving later, before reading it as a different amplitude. That's a suggestion from one radius only, not tested.
