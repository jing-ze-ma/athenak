# tanacc-1107: does the full Eddington tensor (vet_gd_tangential) bring E, F closer to the exact answer?

Branch `tanacc-1107` (fork), commit 36dec1dfb = prod-1107v 698b82413 + the opt-in diagnostic output `m1_lat`.
The output dumps the closure tensor handed to the implicit operator (VlatTT slots, rotate-averaged), the cell unit
vector, and J, H_r of the vet_gd formal solution. Binary: viper GPU build (ROCm 7.2, PROBLEM none),
`bin/athena_tanacc`, md5 5c3cfb3ba14747e1fd75453959d34bc9. All runs: viper apudev, 1 MI300A per run.
No production restart is used: this is a fresh analytic test, so no restart rotation applies.

## 1. Setup

- Problem: built-in `rad_m1_beam`, `m1_test = sph_shell`.
  - Black-body inner sphere r = 1: `implicit_bc_x1min = flux`, F_in = 1, `vet_bc_inner = emit`. In E units with c = 1 this gives I_s = 4.
  - r in [1, 6]. Gas held (`gas_hold`) and cold, so it does not emit. Background chi = 1e-6 (thin).
  - Marshak top with E_bath = 0. Grey, no scattering.
  - A spherical clump of radius a = 0.4 at r = 3 on the axis (sp: theta = pi/2, phi = pi/4; cs: theta = pi/2, phi = 0, at the centre of panel +x).
- Cases (clump extinction chi, source B = arad T^4):

  | case | chi | tau across | B | what it tests |
  |---|---|---|---|---|
  | star | none | - | - | symmetric control (exact anisotropy = 0) |
  | shd | 1000 | 800 | 0 | opaque shadow (umbra/penumbra) |
  | semi | 1.25 | 1 | 0 | tau ~ 1 absorber |
  | mid | 20 | 16 | 0 | thick cold absorber (per-cell tau ~ 1-3) |
  | hot5 | 5 | 4 | 0.1 | emitting, tau ~ few |
  | hot | 50 | 40 | 0.1 | emitting, thick (the interior is in the diffusion limit, E = B) |

  arad = 1e16 with a tiny T_clump keeps the held gas's sound speed from limiting dt.
- Grids (r uniform, 1 MeshBlock along r):
  - sp wedge: theta pi/2 +- 0.6 (periodic + sp_x2_periodic_image), phi [0, pi/2] (periodic). lo 32x36x48, hi 64x72x96.
  - cs: full sphere. lo 32 x 48^2 per panel, hi 64 x 96^2. The clump subtends 0.134 rad: about 4 / 8 cells (lo / hi).
- Radiation keys: the production thin-region recipe (deck generator `mkdecks.py`).
  - Closure and directions: vet_col + vet_gd + vet_gd_replace, taucut 30, geodesic 42 directions (`vet_gd_geo = 2`) with `vet_rotate_avg = 8`, rotate every 2 steps, `vet_col_lat_offdiag = operator`.
  - Face flux: `implicit_flux = blend` idort_f, `implicit_flux_faces = all`, `implicit_flux_beam = halfrange`, `implicit_thin_switch = fse` + `hc`, plus the final_1012 thin keys.
  - Solver: precond mg, time_scheme be, implicit_cfl 1000. 60 cycles = 16-350 light crossings: steady.
  - Only `vet_gd_tangential` differs between arms.
  - The high-direction arm "+hd": sp GL 18x24 = 432 directions (invariant at the wedge walls), cs geodesic level 6 = 362 directions, both x rotate_avg 8.
- Decks: `inp/`; runs: `runs/<g>_<lev>_<case>_<on|off>[_hd]`.
- Smokes: 12804243 (found the cs precond requirement) and 12805594 (clean, the decks as run). mid/hot5 smoke: 12814542.
- Production jobs: 12806889, 12806893, 12806894, 12806895, 12817540, 12817541, 12817542. Every run reached cycle 60, rc=0, no FATAL, no non-convergence.

## 2. Analytic reference (`ana/exact.py`)

I(n) = I_s [the backward ray hits the star] exp(-tau(n)) + B (1 - exp(-tau(n))). tau(n) is chi times the chord through the
clump, clipped at the star. J, H, K come from product Gauss-Legendre x uniform-phi quadrature: 96x160 on the star cap,
and the same on the clump cap (the full sphere for cells inside the clump).

Self-checks:
- Star only: J = 2(1 - sqrt(1 - 1/r^2)), H r^2 = 1 and D_rr = (1 + mu0 + mu0^2)/3, all reproduced to 1e-8.
- Shadow J at the penumbra points changes by 1e-3 relative from 96 to 192 mu nodes.
- Hot clump centre: J = B = 0.1 exactly.

Evaluation cells (`ana/evalrun.py`): the radial shells nearest r = 2.0, 3.6, 4.5, 5.5 within 0.5 rad of the clump axis, plus
every cell inside an emitting clump. Errors are normalised by the unshadowed field J0(r) = 2(1 - sqrt(1-1/r^2)) for E, J_fs and
the pressure anisotropy, and by 1/r^2 for F. Regions: "free" means |J_ex/J0 - 1| < 0.02; "aff" (clump-affected) otherwise;
"aff_far" is aff with r > 4.

Quantities, all against the EXACT solution:
- E_L1 = <|E - J_ex|>/J0.
- Fr_L1 = <|F_r - H_r,ex|> r^2. Flat (sp only) is the lateral flux.
- Q = the trace-free tangential pressure anisotropy ((P_tt - P_pp)/2, P_tp): sp in full components; cs as a magnitude, because the cs frame is cs_qf.
- Q_used_err: the error of the anisotropy the operator actually uses. That is 0 with tan off, so the error equals |Q_ex|.
- Drr_L1: the error of the handed-over D_rr, times J/J0.
- Jfs_L1: the error of the current-set vet_gd formal-solution J (the "SC J").

## 3. Results

### 3.1 On vs off, cell by cell (`ana/ondiff.py`, every cell)

| case | sp lo / hi / lo+hd | cs lo / hi / lo+hd |
|---|---|---|
| star, shd, semi, mid | BITWISE identical E, F in all arms | BITWISE identical in all arms |
| hot5 | max dE/J0 1.5e-2 / 1.5e-2 / -, mean 6e-5 / 8e-5 | 1.9e-2 / 1.7e-2 / -, mean 8e-6 / 1.1e-5 |
| hot | max dE/J0 1.8e-2 / 2.8e-3 / 4.7e-2, mean 1.2e-5 / 6.5e-6 / 2.2e-5 | 6.2e-2 / 3.6e-3 / 5.8e-2, mean 7e-6 / 1e-6 / 4e-6 |

With the production blend + fse/hc thin switch, the face flux in every cell where the field streams (high |H|/J) is the
half-range formal-solution flux, not the tensor flux. This holds even inside an opaque or tau-16 COLD absorber, where
the beam still streams. So the tensor, and with it vet_gd_tangential, enters only in emitting, near-isotropic cells: the
hot clumps and their surfaces. In every thin or absorbing case it is an exact no-op. The cs log reports
"Q active in 20 of 20 transverse assemblies", but the result is bitwise the same.

### 3.2 Where the key does act (emitting clumps), errors vs exact on the cells where on and off differ by > 1e-4 J0 (`ana/diffreg.py`)

Error normalised by max(J0, J_ex). "on closer" = share of those cells where |E_on - J_ex| < |E_off - J_ex|.

| run | region | N | E err on / off | on closer | r^2 F_r err on / off |
|---|---|---|---|---|---|
| sp_lo_hot | inside clump | 84 | 5.763e-02 / 5.760e-02 | 46 % | 1.723e-01 / 1.723e-01 |
| sp_lo_hot | shell 1-1.5 a | 103 | 1.424e-01 / 1.419e-01 | 17 % | 9.168e-02 / 9.151e-02 |
| sp_hi_hot | inside clump | 582 | 6.312e-02 / 6.321e-02 | 56 % | 1.090e-01 / 1.091e-01 |
| sp_hi_hot | shell 1-1.5 a | 996 | 1.426e-01 / 1.423e-01 | 16 % | 1.185e-01 / 1.183e-01 |
| sp_lo_hot_hd | inside clump | 90 | 7.118e-02 / 7.118e-02 | 48 % | 1.728e-01 / 1.728e-01 |
| sp_lo_hot_hd | shell 1-1.5 a | 138 | 1.158e-01 / 1.160e-01 | 19 % | 6.888e-02 / 6.878e-02 |
| cs_lo_hot | inside clump | 101 | 7.813e-02 / 7.808e-02 | 43 % | 1.704e-01 / 1.704e-01 |
| cs_lo_hot | shell 1-1.5 a | 161 | 1.110e-01 / 1.085e-01 | 28 % | 7.599e-02 / 7.556e-02 |
| cs_hi_hot | inside clump | 601 | 6.507e-02 / 6.514e-02 | 49 % | 9.618e-02 / 9.621e-02 |
| cs_hi_hot | shell 1-1.5 a | 1132 | 9.951e-02 / 9.928e-02 | 31 % | 8.849e-02 / 8.832e-02 |
| cs_lo_hot_hd | inside clump | 94 | 6.495e-02 / 6.491e-02 | 40 % | 1.708e-01 / 1.708e-01 |
| cs_lo_hot_hd | shell 1-1.5 a | 153 | 1.190e-01 / 1.177e-01 | 6 % | 7.102e-02 / 7.077e-02 |
| sp_lo_hot5 | inside clump | 172 | 8.858e-02 / 8.756e-02 | 35 % | 6.415e-02 / 6.411e-02 |
| sp_lo_hot5 | shell 1-1.5 a | 303 | 9.062e-02 / 8.986e-02 | 25 % | 7.048e-02 / 6.989e-02 |
| sp_hi_hot5 | inside clump | 1378 | 6.479e-02 / 6.350e-02 | 36 % | 3.809e-02 / 3.847e-02 |
| sp_hi_hot5 | shell 1-1.5 a | 1859 | 8.425e-02 / 8.278e-02 | 28 % | 7.032e-02 / 6.916e-02 |
| cs_lo_hot5 | inside clump | 184 | 1.057e-01 / 1.051e-01 | 36 % | 7.436e-02 / 7.439e-02 |
| cs_lo_hot5 | shell 1-1.5 a | 317 | 8.824e-02 / 8.762e-02 | 40 % | 7.459e-02 / 7.424e-02 |
| cs_hi_hot5 | inside clump | 1408 | 6.948e-02 / 6.834e-02 | 40 % | 4.273e-02 / 4.287e-02 |
| cs_hi_hot5 | shell 1-1.5 a | 1881 | 7.014e-02 / 6.890e-02 | 34 % | 6.464e-02 / 6.384e-02 |

(The full lines, including the "all" and "outside 1.5 a" regions, are in `ana/diffreg_hot.txt` and `ana/diffreg_E.txt`.)

In the emitting arms, sp and cs, lo, hi and +hd, the full tensor mostly moves E and F slightly AWAY from the exact solution.
- The mean change is 0.1-2 % of the error itself.
- On the clump surface shell (1-1.5 a), the anisotropic region, only 6-40 % of the changed cells get closer, in all 10 arms.
- Inside the clumps the result is mixed: 35-56 % of cells get closer, and the means differ by < 2 %.
- The size of the effect is far below the scheme's error: E is 6-12 % off on the clump and its surface either way.

### 3.3 All cases, clump-affected region (star: all cells); on / off (`ana/summary_rows.md`)

| arm | region N | E_L1 | E_bias | Fr_L1 | Q_ex | Q_used_err | Drr_L1 | Jfs_L1 |
|---|---|---|---|---|---|---|---|---|
| sp lo star | all 2912 | 0.1363 / 0.1363 | +0.0719 / +0.0719 | 0.1325 / 0.1325 | 5.3e-10 | 5.6e-03 / 5.3e-10 | 9.6e-03 | 0.451 |
| sp hi star | all 11680 | 0.1566 / 0.1566 | +0.0829 / +0.0829 | 0.1527 / 0.1527 | 5.6e-10 | 6.1e-03 / 5.6e-10 | 9.5e-03 | 0.490 |
| sp lo shd | aff 496 | 0.1483 / 0.1483 | -0.0709 / -0.0709 | 0.1479 / 0.1479 | 1.1e-03 | 3.1e-03 / 1.1e-03 | 8.8e-03 | 0.450 |
| sp hi shd | aff 2004 | 0.1559 / 0.1559 | -0.0541 / -0.0541 | 0.1556 / 0.1556 | 1.1e-03 | 3.3e-03 / 1.1e-03 | 7.9e-03 | 0.526 |
| sp lo+hd shd | aff 496 | 0.0454 / 0.0454 | -0.0328 / -0.0328 | 0.0459 / 0.0459 | 1.1e-03 | 6.8e-04 / 1.1e-03 | 3.1e-03 | 0.060 |
| sp lo semi | aff 420 | 0.1531 / 0.1531 | -0.0583 / -0.0583 | 0.1513 / 0.1513 | 6.2e-04 | 3.6e-03 / 6.2e-04 | 1.0e-02 | 0.507 |
| sp hi semi | aff 1680 | 0.1652 / 0.1652 | -0.0560 / -0.0560 | 0.1643 / 0.1643 | 6.3e-04 | 3.9e-03 / 6.3e-04 | 9.3e-03 | 0.565 |
| sp lo mid | aff 496 | 0.1363 / 0.1363 | -0.0436 / -0.0436 | 0.1352 / 0.1352 | 1.1e-03 | 3.2e-03 / 1.1e-03 | 8.8e-03 | 0.452 |
| sp hi mid | aff 1996 | 0.1428 / 0.1428 | -0.0326 / -0.0326 | 0.1421 / 0.1421 | 1.1e-03 | 3.4e-03 / 1.1e-03 | 8.0e-03 | 0.526 |
| sp lo hot5 | aff 1132 | 0.1174 / 0.1171 | +0.0180 / +0.0181 | 0.1125 / 0.1124 | 4.9e-03 | 5.5e-03 / 4.9e-03 | 1.5e-02 | 0.315 |
| sp hi hot5 | aff 5204 | 0.1198 / 0.1193 | +0.0328 / +0.0329 | 0.1109 / 0.1108 | 4.2e-03 | 5.5e-03 / 4.2e-03 | 1.4e-02 | 0.327 |
| sp lo hot | aff 1344 | 0.1154 / 0.1154 | +0.0232 / +0.0233 | 0.1233 / 0.1233 | 5.2e-03 | 5.3e-03 / 5.2e-03 | 1.2e-02 | 0.329 |
| sp hi hot | aff 5960 | 0.1210 / 0.1210 | +0.0297 / +0.0298 | 0.1246 / 0.1246 | 4.7e-03 | 4.6e-03 / 4.7e-03 | 1.1e-02 | 0.342 |
| sp lo+hd hot | aff 1344 | 0.0307 / 0.0307 | -0.0091 / -0.0091 | 0.0380 / 0.0380 | 5.2e-03 | 2.1e-03 / 5.2e-03 | 7.4e-03 | 0.052 |
| cs lo star | all 3056 | 0.1822 / 0.1822 | -0.0056 / -0.0056 | 0.1805 / 0.1805 | 4.7e-10 | 5.2e-03 / 4.7e-10 | 6.5e-03 | 0.510 |
| cs hi star | all 12176 | 0.1966 / 0.1966 | -0.0044 / -0.0044 | 0.1957 / 0.1957 | 5.1e-10 | 5.3e-03 / 5.1e-10 | 5.8e-03 | 0.559 |
| cs lo shd | aff 508 | 0.2512 / 0.2512 | -0.0401 / -0.0401 | 0.2523 / 0.2523 | 1.1e-03 | 3.0e-03 / 1.1e-03 | 5.7e-03 | 0.515 |
| cs hi shd | aff 2064 | 0.2646 / 0.2646 | -0.0288 / -0.0288 | 0.2656 / 0.2656 | 1.1e-03 | 3.2e-03 / 1.1e-03 | 4.9e-03 | 0.588 |
| cs lo+hd shd | aff 508 | 0.0475 / 0.0475 | -0.0331 / -0.0331 | 0.0482 / 0.0482 | 1.1e-03 | 4.7e-04 / 1.1e-03 | 2.8e-03 | 0.050 |
| cs lo semi | aff 424 | 0.2580 / 0.2580 | -0.0075 / -0.0075 | 0.2585 / 0.2585 | 6.2e-04 | 4.2e-03 / 6.2e-04 | 7.0e-03 | 0.584 |
| cs hi semi | aff 1724 | 0.2657 / 0.2657 | -0.0121 / -0.0121 | 0.2665 / 0.2665 | 6.3e-04 | 4.2e-03 / 6.3e-04 | 6.0e-03 | 0.637 |
| cs lo mid | aff 508 | 0.2475 / 0.2475 | -0.0141 / -0.0141 | 0.2487 / 0.2487 | 1.1e-03 | 3.3e-03 / 1.1e-03 | 5.7e-03 | 0.516 |
| cs hi mid | aff 2060 | 0.2624 / 0.2624 | -0.0074 / -0.0074 | 0.2630 / 0.2630 | 1.1e-03 | 3.4e-03 / 1.1e-03 | 5.0e-03 | 0.589 |
| cs lo hot5 | aff 1172 | 0.1662 / 0.1659 | -0.0076 / -0.0075 | 0.1611 / 0.1609 | 4.9e-03 | 3.1e-03 / 4.9e-03 | 1.1e-02 | 0.363 |
| cs hi hot5 | aff 5436 | 0.1618 / 0.1614 | +0.0088 / +0.0089 | 0.1550 / 0.1549 | 4.2e-03 | 3.1e-03 / 4.2e-03 | 9.8e-03 | 0.369 |
| cs lo hot | aff 1364 | 0.1607 / 0.1605 | -0.0174 / -0.0171 | 0.1694 / 0.1692 | 5.3e-03 | 3.1e-03 / 5.3e-03 | 9.0e-03 | 0.367 |
| cs hi hot | aff 6156 | 0.1614 / 0.1614 | -0.0103 / -0.0102 | 0.1663 / 0.1663 | 4.7e-03 | 2.9e-03 / 4.7e-03 | 9.0e-03 | 0.379 |
| cs lo+hd hot | aff 1364 | 0.0303 / 0.0302 | -0.0186 / -0.0185 | 0.0370 / 0.0369 | 5.3e-03 | 1.9e-03 / 5.3e-03 | 5.4e-03 | 0.036 |

Columns with a single value are identical in both arms (the tensor slots and the formal solution do not depend on the key).

### 3.4 Accuracy of the handed-over tensor itself (exact vs code; Q_used_err column)

- The exact tangential anisotropy is small everywhere in these tests.
  - Q_ex / J0 = 1e-3 in shadow, semi-transparent and mid regions.
  - 4-5e-3 around the hot clumps.
  - 0 for the star.
- With production directions (42 x rotate_avg 8), the ray noise in the handed-over Q is 3-6e-3 J0. That is 3-5 times the signal in the absorber cases. Turning the key on raises the tensor error there: Q_used_err 3.1-3.6e-3 on vs 1.0-1.1e-3 off.
  - Near the hot clumps on and off are about equal on sp (aff: 4.6e-3 / 4.7e-3 hi).
  - cs is better on (2.9e-3 vs 4.7e-3).
- With the high-direction sets, the handed-over Q is more accurate than zero in the affected regions:
  - sp shd 6.8e-4 vs 1.1e-3; sp hot 2.1e-3 vs 5.2e-3;
  - cs shd 4.7e-4 vs 1.1e-3; cs hot 1.9e-3 vs 5.3e-3.
- So the tensor component is right only when the direction set resolves it. Even then, E and F do not improve (3.1, 3.2).
- D_rr (key-independent): error 0.6-1.5 % (geo2), 0.3-0.7 % (+hd).
- D_r,lat (key-independent): `vlat_err` in `ana/table_blend.md`, 3-5e-2 J0 with geo2 and 4-8e-3 with +hd.

### 3.5 Convergence

- Spatial (lo -> hi, factor 2 in every direction): no convergence with the production direction set.
  - E_L1 grows slightly: sp star 0.136 -> 0.157, cs star 0.182 -> 0.197, sp shd aff 0.148 -> 0.156.
  - This is the known effect that more shells mean more lateral interpolations of a ray-sampled field (sptests T6a).
- Directions (42 x 8 -> 432 or 362 x 8, at lo): the error falls about 10x.
  - sp shd free 0.125 -> 0.012, cs shd free 0.162 -> 0.011; aff 0.148 -> 0.045 (sp), 0.251 -> 0.048 (cs).
  - Fr follows the same pattern.
- So the error of E and F in thin, streaming regions is set by the direction set feeding the half-range flux, not by the closure tensor.
- The key's effect does not change with either refinement: bitwise no-op in non-emitting cases, slightly adverse in emitting ones.

### 3.6 Comparison with the code's own formal solution (the "SC J")

The current-set J_fs is 30-70 % off the exact J at the sampled shells with 42 directions (column Jfs_L1). It is 3-6 % off
with the high-direction sets. The moment E, which uses the rotate-averaged closure, is 3-4x closer to the exact answer than
the single-set J_fs.

So "E vs SC J" (Caltech cstan_1106: E 2-4 % further from SC J with the tensor on) is not an accuracy measure in thin regions
unless SC J is itself converged in angle. Against the exact solution here, the direction of the change agrees with Caltech:
on is slightly worse. The size is 10-100x smaller than 2-4 %, because these test states have only small emitting regions.

### 3.7 Inconclusive / not done

- Pure-VET arm (`implicit_flux = central`, decks `inp/*_C`, runs `runsC/`): this arm would make the tensor carry the thin transport.
  - The inner flux boundary injects nothing on that path: E ~ 5e-11 everywhere, J_fs = 0.
  - So it tests nothing. It is not a production path. Noted as a possible separate bug or limitation.
- No case with a strongly anisotropic, near-isotropic emitting region larger than a few cells, i.e. a resolved tau ~ 1 photosphere with lateral structure.
  - The hot clumps are the only regime where the key acts, and they are 4-8 cells across.
  - A stellar photosphere with granulation is the production case. It has no analytic reference.
- Seams: the cs clump sits at a panel centre. Seam effects were not tested.
- The sp wedge is periodic in theta: cells within ~0.3 rad of the theta walls are excluded (only cells within 0.5 rad of the axis are evaluated).

## 4. Verdict

No. On both sp and cs, the full tensor does not bring E or F closer to the exact solution in any tested regime.
- Thin, shadowed, tau~1 and thick cold absorbers: with the production blend/half-range flux, E and F are bitwise identical with the key on or off.
- Emitting clumps (tau 4 and 40): the key changes E and F by up to 2-6 % of J0 locally, 1e-5 on average. The change moves them slightly AWAY from the exact answer: the error grows by 0.1-2 %, and 60-94 % of the changed surface cells get worse. Inside the clumps the result is mixed (35-56 % of cells closer).
- The tensor components are more accurate only with a high direction count. With the production 42 x 8 set, the handed-over tangential anisotropy is mostly ray noise (3-5x the exact signal).
- The dominant error in these tests is the direction set: 13-30 % with 42 x 8, 1-5 % with ~400 x 8. It is not the closure tensor.

## 5. Files

- Decks / generator: `mkdecks.py`, `inp/`. Job script: `jobs/runlist.sbatch`. Logs: `logs/`.
- Analysis: `ana/exact.py`, `ana/evalrun.py` (one `runs/*/eval.json` per run), `ana/table.py` (`ana/table_blend.md` = every metric and region), `ana/ondiff.py` (`ana/ondiff.txt`), `ana/diffreg.py` (`ana/diffreg_hot.txt`, `ana/diffreg_E.txt`), `ana/summary_rows.md`.
